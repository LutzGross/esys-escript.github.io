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

#include <sc_private.h>
#include <sc_refcount.h>

void
sc_refcount_init_invalid (sc_refcount_t * rc)
{
  SC_ASSERT (rc != NULL);

  rc->package_id = -1;
  rc->refcount = -1;
}

void
sc_refcount_init (sc_refcount_t * rc, int package_id)
{
  SC_ASSERT (rc != NULL);
  SC_ASSERT (package_id == -1 || sc_package_is_registered (package_id));

  rc->package_id = package_id;
  rc->refcount = 1;

#ifdef SC_ENABLE_DEBUG
  sc_package_rc_count_add (rc->package_id, 1);
#endif
}

sc_refcount_t      *
sc_refcount_new (int package_id)
{
  sc_refcount_t      *rc;

  rc = SC_ALLOC (sc_refcount_t, 1);
  sc_refcount_init (rc, package_id);

  return rc;
}

void
sc_refcount_destroy (sc_refcount_t * rc)
{
  SC_ASSERT (rc != NULL);
  SC_ASSERT (rc->refcount == 0);

  SC_FREE (rc);
}

void
sc_refcount_ref (sc_refcount_t * rc)
{
  SC_ASSERT (rc != NULL);
  SC_ASSERT (rc->refcount > 0);

  ++rc->refcount;
}

int
sc_refcount_unref (sc_refcount_t * rc)
{
  SC_ASSERT (rc != NULL);
  SC_ASSERT (rc->refcount > 0);

  if (--rc->refcount == 0) {
#ifdef SC_ENABLE_DEBUG
    sc_package_rc_count_add (rc->package_id, -1);
#endif
    return 1;
  }
  else {
    return 0;
  }
}

int
sc_refcount_is_active (const sc_refcount_t * rc)
{
  SC_ASSERT (rc != NULL);

  return rc->refcount > 0;
}

int
sc_refcount_is_last (const sc_refcount_t * rc)
{
  SC_ASSERT (rc != NULL);

  return rc->refcount == 1;
}
