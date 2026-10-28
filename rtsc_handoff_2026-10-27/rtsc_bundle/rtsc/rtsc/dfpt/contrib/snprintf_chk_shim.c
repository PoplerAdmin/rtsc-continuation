/* Workaround for Ubuntu 24.04 quantum-espresso 6.7 (Debian build with _FORTIFY_SOURCE):
 * clib/md5_from_file.c calls snprintf(&buf[2*i], 8, "%02x", b) on a 33-byte buffer, so
 * near the end of the buffer the declared maxlen (8) exceeds the remaining object size and
 * glibc's __snprintf_chk aborts ("buffer overflow detected") although only 3 bytes are written.
 * This shim clamps maxlen to the object size (never writes beyond what fortify allows) instead
 * of aborting. It does not change any numerical code path. */
#include <stdarg.h>
#include <stdio.h>
#include <stddef.h>
int __snprintf_chk(char *s, size_t maxlen, int flag, size_t slen, const char *fmt, ...) {
    (void)flag;
    if (maxlen > slen) maxlen = slen;
    va_list ap; va_start(ap, fmt);
    int r = vsnprintf(s, maxlen, fmt, ap);
    va_end(ap);
    return r;
}
