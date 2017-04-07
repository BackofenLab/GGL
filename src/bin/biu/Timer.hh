#ifndef TIMER_HH_
#define TIMER_HH_

#include <ctime>

namespace biu
{

		/**
		 * Timer class to measure runtime in milliseconds.
		 *
		 * @author Martin Mann <mmann@@informatik.uni-freiburg.de>
		 */
	class Timer {
		private:
				//! starting time
			clock_t t0;
		public:
				//! Sets starting time.
			void start(void){
				t0 = clock();
			}
				//! Returns time consumption in milliseconds until now from last
				//! start() call on.
			double stop(void) {
				return (static_cast<double>(clock()-t0) / CLOCKS_PER_SEC) * 1000.0;
			}
	};
	
} // namespace biu
#endif /*TIMER_HH_*/
