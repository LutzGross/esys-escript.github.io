/*
  This file is part of the SC Library.
  The SC Library provides support for parallel scientific applications.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors

  The SC Library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  The SC Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with the SC Library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
*/

#include <sc_io.h>
#include <libb64.h>
#ifdef SC_HAVE_ZLIB
#include <zlib.h>
#endif

sc_io_sink_t       *
sc_io_sink_new (sc_io_type_t iotype, sc_io_mode_t mode,
                sc_io_encode_t encode, ...)
{
  sc_io_sink_t       *sink;
  va_list             ap;

  SC_ASSERT (0 <= iotype && iotype < SC_IO_TYPE_LAST);
  SC_ASSERT (0 <= mode && mode < SC_IO_MODE_LAST);
  SC_ASSERT (0 <= encode && encode < SC_IO_ENCODE_LAST);

  sink = SC_ALLOC_ZERO (sc_io_sink_t, 1);
  sink->iotype = iotype;
  sink->mode = mode;
  sink->encode = encode;

  va_start (ap, encode);
  if (iotype == SC_IO_TYPE_BUFFER) {
    sink->buffer = va_arg (ap, sc_array_t *);
    if (sink->mode == SC_IO_MODE_WRITE) {
      sc_array_resize (sink->buffer, 0);
    }
  }
  else if (iotype == SC_IO_TYPE_FILENAME) {
    const char         *filename = va_arg (ap, const char *);

    sink->file = fopen (filename,
                        sink->mode == SC_IO_MODE_WRITE ? "wb" : "ab");
    if (sink->file == NULL) {
      SC_FREE (sink);
      return NULL;
    }
  }
  else if (iotype == SC_IO_TYPE_FILEFILE) {
    sink->file = va_arg (ap, FILE *);
    if (ferror (sink->file)) {
      SC_FREE (sink);
      return NULL;
    }
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
  va_end (ap);

  return sink;
}

int
sc_io_sink_destroy (sc_io_sink_t * sink)
{
  int                 retval;

  /* The error value SC_IO_ERROR_AGAIN is turned into FATAL */
  retval = sc_io_sink_complete (sink, NULL, NULL);
  if (sink->iotype == SC_IO_TYPE_FILENAME) {
    SC_ASSERT (sink->file != NULL);

    /* Attempt close even on complete error */
    retval = fclose (sink->file) || retval;
  }
  SC_FREE (sink);

  return retval ? SC_IO_ERROR_FATAL : SC_IO_ERROR_NONE;
}

int
sc_io_sink_write (sc_io_sink_t * sink, const void *data, size_t bytes_avail)
{
  size_t              bytes_out;

  bytes_out = 0;

  if (sink->iotype == SC_IO_TYPE_BUFFER) {
    size_t              elem_size, new_count;

    SC_ASSERT (sink->buffer != NULL);
    elem_size = sink->buffer->elem_size;
    new_count =
      (sink->buffer_bytes + bytes_avail + elem_size - 1) / elem_size;
    sc_array_resize (sink->buffer, new_count);
    /* For a view sufficient size is asserted only in debug mode. */
    if (new_count * elem_size > SC_ARRAY_BYTE_ALLOC (sink->buffer)) {
      return SC_IO_ERROR_FATAL;
    }

    memcpy (sink->buffer->array + sink->buffer_bytes, data, bytes_avail);
    sink->buffer_bytes += bytes_avail;
    bytes_out = bytes_avail;
  }
  else if (sink->iotype == SC_IO_TYPE_FILENAME ||
           sink->iotype == SC_IO_TYPE_FILEFILE) {
    SC_ASSERT (sink->file != NULL);
    bytes_out = fwrite (data, 1, bytes_avail, sink->file);
    if (bytes_out != bytes_avail) {
      return SC_IO_ERROR_FATAL;
    }
  }

  sink->bytes_in += bytes_avail;
  sink->bytes_out += bytes_out;

  return SC_IO_ERROR_NONE;
}

int
sc_io_sink_complete (sc_io_sink_t * sink,
                     size_t * bytes_in, size_t * bytes_out)
{
  int                 retval;

  retval = 0;
  if (sink->iotype == SC_IO_TYPE_BUFFER) {
    SC_ASSERT (sink->buffer != NULL);
    if (sink->buffer_bytes % sink->buffer->elem_size != 0) {
      return SC_IO_ERROR_AGAIN;
    }
  }
  else if (sink->iotype == SC_IO_TYPE_FILENAME ||
           sink->iotype == SC_IO_TYPE_FILEFILE) {
    SC_ASSERT (sink->file != NULL);
    retval = fflush (sink->file);
  }
  if (retval) {
    return SC_IO_ERROR_FATAL;
  }

  if (bytes_in != NULL) {
    *bytes_in = sink->bytes_in;
  }
  if (bytes_out != NULL) {
    *bytes_out = sink->bytes_out;
  }
  sink->bytes_in = sink->bytes_out = 0;

  return SC_IO_ERROR_NONE;
}

int
sc_io_sink_align (sc_io_sink_t * sink, size_t bytes_align)
{
  size_t              fill_bytes;
  char               *fill;
  int                 retval;

  fill_bytes = (bytes_align - sink->bytes_out % bytes_align) % bytes_align;
  fill = SC_ALLOC_ZERO (char, fill_bytes);
  retval = sc_io_sink_write (sink, fill, fill_bytes);
  SC_FREE (fill);

  return retval;
}

sc_io_source_t     *
sc_io_source_new (sc_io_type_t iotype, sc_io_encode_t encode, ...)
{
  sc_io_source_t     *source;
  va_list             ap;

  SC_ASSERT (0 <= iotype && iotype < SC_IO_TYPE_LAST);
  SC_ASSERT (0 <= encode && encode < SC_IO_ENCODE_LAST);

  source = SC_ALLOC_ZERO (sc_io_source_t, 1);
  source->iotype = iotype;
  source->encode = encode;

  va_start (ap, encode);
  if (iotype == SC_IO_TYPE_BUFFER) {
    source->buffer = va_arg (ap, sc_array_t *);
  }
  else if (iotype == SC_IO_TYPE_FILENAME) {
    const char         *filename = va_arg (ap, const char *);

    source->file = fopen (filename, "rb");
    if (source->file == NULL) {
      SC_FREE (source);
      return NULL;
    }
  }
  else if (iotype == SC_IO_TYPE_FILEFILE) {
    source->file = va_arg (ap, FILE *);
    if (ferror (source->file)) {
      SC_FREE (source);
      return NULL;
    }
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
  va_end (ap);

  return source;
}

int
sc_io_source_destroy (sc_io_source_t * source)
{
  int                 retval;

  /* complete reading */
  retval = sc_io_source_complete (source, NULL, NULL);

  /* destroy mirror */
  if (source->mirror != NULL) {
    retval = sc_io_sink_destroy (source->mirror) || retval;
    sc_array_destroy (source->mirror_buffer);
  }

  /* The error value SC_IO_ERROR_AGAIN is turned into FATAL */
  if (source->iotype == SC_IO_TYPE_FILENAME) {
    SC_ASSERT (source->file != NULL);

    /* Attempt close even on complete error */
    retval = fclose (source->file) || retval;
  }
  SC_FREE (source);

  return retval ? SC_IO_ERROR_FATAL : SC_IO_ERROR_NONE;
}

int
sc_io_source_read (sc_io_source_t * source, void *data,
                   size_t bytes_avail, size_t * bytes_out)
{
  int                 retval;
  size_t              bbytes_out;

  retval = 0;
  bbytes_out = 0;

  if (source->iotype == SC_IO_TYPE_BUFFER) {
    SC_ASSERT (source->buffer != NULL);
    bbytes_out = SC_ARRAY_BYTE_ALLOC (source->buffer);
    SC_ASSERT (bbytes_out >= source->buffer_bytes);
    bbytes_out -= source->buffer_bytes;
    bbytes_out = SC_MIN (bbytes_out, bytes_avail);

    if (data != NULL) {
      memcpy (data, source->buffer->array + source->buffer_bytes, bbytes_out);
    }
    source->buffer_bytes += bbytes_out;
  }
  else if (source->iotype == SC_IO_TYPE_FILENAME ||
           source->iotype == SC_IO_TYPE_FILEFILE) {
    SC_ASSERT (source->file != NULL);
    if (data != NULL) {
      bbytes_out = fread (data, 1, bytes_avail, source->file);
      if (bbytes_out < bytes_avail) {
        retval = !feof (source->file) || ferror (source->file);
      }
      if (retval == SC_IO_ERROR_NONE && source->mirror != NULL) {
        retval = sc_io_sink_write (source->mirror, data, bbytes_out);
      }
    }
    else {
      retval = fseek (source->file, (long) bytes_avail, SEEK_CUR);
      bbytes_out = bytes_avail;
    }
  }
  if (retval) {
    return SC_IO_ERROR_FATAL;
  }
  if (bytes_out == NULL && bbytes_out < bytes_avail) {
    return SC_IO_ERROR_FATAL;
  }

  if (bytes_out != NULL) {
    *bytes_out = bbytes_out;
  }
  source->bytes_in += bbytes_out;
  source->bytes_out += bbytes_out;

  return SC_IO_ERROR_NONE;
}

int
sc_io_source_complete (sc_io_source_t * source,
                       size_t * bytes_in, size_t * bytes_out)
{
  int                 retval = SC_IO_ERROR_NONE;

  if (source->iotype == SC_IO_TYPE_BUFFER) {
    SC_ASSERT (source->buffer != NULL);
    if (source->buffer_bytes % source->buffer->elem_size != 0) {
      return SC_IO_ERROR_AGAIN;
    }
  }
  else if (source->iotype == SC_IO_TYPE_FILENAME ||
           source->iotype == SC_IO_TYPE_FILEFILE) {
    if (source->mirror != NULL) {
      retval = sc_io_sink_complete (source->mirror, NULL, NULL);
    }
  }

  if (bytes_in != NULL) {
    *bytes_in = source->bytes_in;
  }
  if (bytes_out != NULL) {
    *bytes_out = source->bytes_out;
  }
  source->bytes_in = source->bytes_out = 0;

  return retval;
}

int
sc_io_source_align (sc_io_source_t * source, size_t bytes_align)
{
  size_t              fill_bytes;

  fill_bytes = (bytes_align - source->bytes_out % bytes_align) % bytes_align;

  return sc_io_source_read (source, NULL, fill_bytes, NULL);
}

int
sc_io_source_activate_mirror (sc_io_source_t * source)
{
  if (source->iotype == SC_IO_TYPE_BUFFER) {
    return SC_IO_ERROR_FATAL;
  }
  if (source->mirror != NULL) {
    return SC_IO_ERROR_FATAL;
  }

  source->mirror_buffer = sc_array_new (sizeof (char));
  source->mirror = sc_io_sink_new (SC_IO_TYPE_BUFFER, SC_IO_MODE_WRITE,
                                   SC_IO_ENCODE_NONE, source->mirror_buffer);

  return (source->mirror != NULL ? SC_IO_ERROR_NONE : SC_IO_ERROR_FATAL);
}

int
sc_io_source_read_mirror (sc_io_source_t * source, void *data,
                          size_t bytes_avail, size_t * bytes_out)
{
  sc_io_source_t     *mirror_src;
  int                 retval;

  if (source->mirror_buffer == NULL) {
    return SC_IO_ERROR_FATAL;
  }

  mirror_src = sc_io_source_new (SC_IO_TYPE_BUFFER, SC_IO_ENCODE_NONE,
                                 source->mirror_buffer);
  retval = (mirror_src != NULL ? SC_IO_ERROR_NONE : SC_IO_ERROR_FATAL);
  retval = retval || sc_io_source_read (mirror_src, data, bytes_avail,
                                        bytes_out);
  if (mirror_src != NULL) {
    retval = sc_io_source_destroy (mirror_src) || retval;
  }

  return retval;
}

int
sc_vtk_write_binary (FILE * vtkfile, char *numeric_data, size_t byte_length)
{
  size_t              chunks, chunksize, remaining, writenow;
  size_t              code_length, base_length;
  uint32_t            int_header;
  char               *base_data;
  base64_encodestate  encode_state;

  /* VTK format used 32bit header info */
  SC_ASSERT (byte_length <= (size_t) UINT32_MAX);

  /* This value may be changed although this is not tested with VTK */
  chunksize = (size_t) 1 << 15; /* 32768 */
  int_header = (uint32_t) byte_length;

  /* Allocate sufficient memory for base64 encoder */
  code_length = 2 * SC_MAX (chunksize, sizeof (int_header));
  code_length = SC_MAX (code_length, 4) + 1;
  base_data = SC_ALLOC (char, code_length);

  base64_init_encodestate (&encode_state);
  base_length =
    base64_encode_block ((char *) &int_header, sizeof (int_header), base_data,
                         &encode_state);
  SC_ASSERT (base_length < code_length);
  base_data[base_length] = '\0';
  (void) fwrite (base_data, 1, base_length, vtkfile);

  chunks = 0;
  remaining = byte_length;
  while (remaining > 0) {
    writenow = SC_MIN (remaining, chunksize);
    base_length = base64_encode_block (numeric_data + chunks * chunksize,
                                       writenow, base_data, &encode_state);
    SC_ASSERT (base_length < code_length);
    base_data[base_length] = '\0';
    (void) fwrite (base_data, 1, base_length, vtkfile);
    remaining -= writenow;
    ++chunks;
  }

  base_length = base64_encode_blockend (base_data, &encode_state);
  SC_ASSERT (base_length < code_length);
  base_data[base_length] = '\0';
  (void) fwrite (base_data, 1, base_length, vtkfile);

  SC_FREE (base_data);
  if (ferror (vtkfile)) {
    return -1;
  }
  return 0;
}

int
sc_vtk_write_compressed (FILE * vtkfile, char *numeric_data,
                         size_t byte_length)
{
#ifdef SC_HAVE_ZLIB
  int                 retval, fseek1, fseek2;
  size_t              iz;
  size_t              blocksize, lastsize;
  size_t              theblock, numregularblocks, numfullblocks;
  size_t              header_entries, header_size;
  size_t              code_length, base_length;
  long                header_pos, final_pos;
  char               *comp_data, *base_data;
  uint32_t           *compression_header;
  uLongf              comp_length;
  base64_encodestate  encode_state;

  /* compute block sizes */
  blocksize = (size_t) (1 << 15);       /* 32768 */
  lastsize = byte_length % blocksize;
  numregularblocks = byte_length / blocksize;
  numfullblocks = numregularblocks + (lastsize > 0 ? 1 : 0);
  header_entries = 3 + numfullblocks;
  header_size = header_entries * sizeof (uint32_t);

  /* allocate compression and base64 arrays */
  code_length = 2 * SC_MAX (blocksize, header_size) + 4 + 1;
  comp_data = SC_ALLOC (char, code_length);
  base_data = SC_ALLOC (char, code_length);

  /* figure out the size of the header and write a dummy */
  compression_header = SC_ALLOC (uint32_t, header_entries);
  compression_header[0] = (uint32_t) numfullblocks;
  compression_header[1] = (uint32_t) blocksize;
  compression_header[2] = (uint32_t)
    (lastsize > 0 || byte_length == 0 ? lastsize : blocksize);
  for (iz = 3; iz < header_entries; ++iz) {
    compression_header[iz] = 0;
  }
  base64_init_encodestate (&encode_state);
  base_length = base64_encode_block ((char *) compression_header,
                                     header_size, base_data, &encode_state);
  base_length +=
    base64_encode_blockend (base_data + base_length, &encode_state);
  SC_ASSERT (base_length < code_length);
  base_data[base_length] = '\0';
  header_pos = ftell (vtkfile);
  (void) fwrite (base_data, 1, base_length, vtkfile);

  /* write the regular data blocks */
  base64_init_encodestate (&encode_state);
  for (theblock = 0; theblock < numregularblocks; ++theblock) {
    comp_length = code_length;
    retval = compress2 ((Bytef *) comp_data, &comp_length,
                        (const Bytef *) (numeric_data + theblock * blocksize),
                        (uLong) blocksize, Z_BEST_COMPRESSION);
    SC_CHECK_ZLIB (retval);
    compression_header[3 + theblock] = comp_length;
    base_length = base64_encode_block (comp_data, comp_length,
                                       base_data, &encode_state);
    SC_ASSERT (base_length < code_length);
    base_data[base_length] = '\0';
    (void) fwrite (base_data, 1, base_length, vtkfile);
  }

  /* write odd-sized last block if necessary */
  if (lastsize > 0) {
    comp_length = code_length;
    retval = compress2 ((Bytef *) comp_data, &comp_length,
                        (const Bytef *) (numeric_data + theblock * blocksize),
                        (uLong) lastsize, Z_BEST_COMPRESSION);
    SC_CHECK_ZLIB (retval);
    compression_header[3 + theblock] = comp_length;
    base_length = base64_encode_block (comp_data, comp_length,
                                       base_data, &encode_state);
    SC_ASSERT (base_length < code_length);
    base_data[base_length] = '\0';
    (void) fwrite (base_data, 1, base_length, vtkfile);
  }

  /* write base64 end block */
  base_length = base64_encode_blockend (base_data, &encode_state);
  SC_ASSERT (base_length < code_length);
  base_data[base_length] = '\0';
  (void) fwrite (base_data, 1, base_length, vtkfile);

  /* seek back, write header block, seek forward */
  final_pos = ftell (vtkfile);
  base64_init_encodestate (&encode_state);
  base_length = base64_encode_block ((char *) compression_header,
                                     header_size, base_data, &encode_state);
  base_length +=
    base64_encode_blockend (base_data + base_length, &encode_state);
  SC_ASSERT (base_length < code_length);
  base_data[base_length] = '\0';
  fseek1 = fseek (vtkfile, header_pos, SEEK_SET);
  (void) fwrite (base_data, 1, base_length, vtkfile);
  fseek2 = fseek (vtkfile, final_pos, SEEK_SET);

  /* clean up and return */
  SC_FREE (compression_header);
  SC_FREE (comp_data);
  SC_FREE (base_data);
  if (fseek1 != 0 || fseek2 != 0 || ferror (vtkfile)) {
    return -1;
  }
#else
  SC_ABORT ("Configure did not find a recent enough zlib.  Abort.\n");
#endif

  return 0;
}

void
sc_fwrite (const void *ptr, size_t size, size_t nmemb, FILE * file,
           const char *errmsg)
{
  size_t              nwritten;

  nwritten = fwrite (ptr, size, nmemb, file);
  SC_CHECK_ABORT (nwritten == nmemb, errmsg);
}

void
sc_fread (void *ptr, size_t size, size_t nmemb, FILE * file,
          const char *errmsg)
{
  size_t              nread;

  nread = fread (ptr, size, nmemb, file);
  SC_CHECK_ABORT (nread == nmemb, errmsg);
}

#ifdef SC_ENABLE_MPIIO

void
sc_mpi_write (MPI_File mpifile, const void *ptr, size_t zcount,
              sc_MPI_Datatype t, const char *errmsg)
{
#ifdef SC_ENABLE_DEBUG
  int                 icount;
#endif
  int                 mpiret;
  sc_MPI_Status       mpistatus;

  mpiret = MPI_File_write (mpifile, (void *) ptr,
                           (int) zcount, t, &mpistatus);
  SC_CHECK_ABORT (mpiret == sc_MPI_SUCCESS, errmsg);

#ifdef SC_ENABLE_DEBUG
  sc_MPI_Get_count (&mpistatus, t, &icount);
  SC_CHECK_ABORT (icount == (int) zcount, errmsg);
#endif
}

#endif
