/* source: http://www.jera.com/techinfo/jtns/jtn002.html */

#ifndef MIN_UNIT_H
#define MIN_UNIT_H

#define mu_assert(message, test) do { if (!(test)) return message; } while (0)
#define mu_run(fn)   do { char *message = fn(); tests_run++; \
                               if (message) return message; } while (0)                                                

extern int tests_run;

/* Add error code checking */  
#define BUFSIZE 64                
#define mu_check(txt,err) do { if(err) { sprintf(buf, "%s %d", txt, err); \
                               return buf; } } while (0)

extern char buf[BUFSIZE];

#endif 