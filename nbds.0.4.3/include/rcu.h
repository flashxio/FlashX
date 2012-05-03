/* 
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 */
#ifndef RCU_H
#define RCU_H

void rcu_update (void);
void rcu_defer_free (void *x);

#endif//RCU_H
