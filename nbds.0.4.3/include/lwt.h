/* 
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 *
 * lightweight tracing
 */ 
#ifndef LWT_H
#define LWT_H

#ifndef ENABLE_TRACE
#define TRACE(...) do { } while (0)
#else
#define TRACE(flag, format, v1, v2) lwt_trace(flag, format, (size_t)(v1), (size_t)(v2))
#endif

#ifndef NDEBUG
#define ASSERT(x) do { if (!(x)) { lwt_halt(); assert(!#x); } } while (0)
#else
#define ASSERT(x) do { } while (0)
#endif

// Dump trace records to <file_name>. The file should be post-processed with "sort" before viewing.
void lwt_dump (const char *file_name) __attribute__ ((externally_visible));

// <flags> indicates what kind of trace messages should be included in the dump. <flags> is a sequence of letters
// followed by numbers (e.g. "x1c9n2g3"). The letters indicate trace categories and the numbers are trace levels 
// for each category. If a category appears in <flags>, then messages from that category will be included in the
// dump if they have a trace level less than or equal to the one specified in <flags>. Categories are case
// sensitive.
void lwt_set_trace_level (const char *flags);

// <flag> is a two character string containing a letter followed by a number (e.g. "f3"). The letter indicates a
// trace category, and the number a trace level. <flag> controls whether or not the trace message gets included in
// the dump. It is only included when its specified category is enabled at a trace level greater than or equal to
// the one in <flag>. Categories are case sensitive. 
static inline void lwt_trace (const char *flag, const char *format, size_t value1, size_t value2) {
    extern char flag_state_[256];
    if (EXPECT_FALSE(flag_state_[(unsigned)flag[0]] >= flag[1])) {
        // embed <flags> in <format> so we don't have to make the lwt_record_t any bigger than it already is
        uint64_t f = ((uint64_t)(size_t)format | ((uint64_t)flag[0] << 56) | ((uint64_t)flag[1] << 48));
        extern void lwt_trace_i (uint64_t format, size_t value1, size_t value2);
        lwt_trace_i(f, value1, value2);
    }
}

void lwt_halt (void);

#endif//LWT_H
