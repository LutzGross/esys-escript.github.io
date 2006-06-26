/* qsortG - drop in replacement for qsort */
/* Useful for consistency across platforms as quicksort isn't guaranteed to be stable (and isn't on Win32) */

#if !defined qsortG_H
#define qsortG_H

extern void qsortG(void *base, size_t nmemb, size_t size,
           int (*compare)(const void *, const void *));

#endif