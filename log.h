#if (defined simple_timer || defined vampir)
#define START_LOG(a) call start_log(a ##_begin, #a )
#define	END_LOG(a) call end_log(a ## _end, #a )
#else
#define START_LOG(a)	! no logging
#define END_LOG(a)	! no logging
#endif

