// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_SYSTEM_EVENT_H
#define SEQAN_HEADER_SYSTEM_EVENT_H

//#include <iterator>

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES EventDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

    struct Event    // this class mustn't exceed the size of HANDLE (needed by waitForAll/Any)
    {
        typedef HANDLE Handle;
		enum { Infinite = INFINITE };
        Handle hEvent;

        Event():
            hEvent(NULL) {}

        Event(BOOL initial) {
            SEQAN_DO_SYS2(open(initial), "Could not create Event");
        }

        ~Event() {
            if (*this) SEQAN_DO_SYS2(close(), "Could not destroy Event");
        }

        Event(Event const &origin) {
            // resource sharing is not yet supported (performance reason)
            // it needs a reference counting technique
            if (origin) {
                hEvent = origin.hEvent;
                const_cast<Event&>(origin).hEvent = NULL;
            } else
                hEvent = NULL;
        }

        inline Event& operator=(Event const &origin) {
            // resource sharing is not yet supported (performance reason)
            // it needs a reference counting technique
            if (origin) {
                hEvent = origin.hEvent;
                const_cast<Event&>(origin).hEvent = NULL;
            } else
                hEvent = NULL;
            return *this;
        }

        inline bool open(BOOL initial = FALSE) {
            return (hEvent = CreateEvent(&EventDefaultAttributes, TRUE, initial, NULL)) != NULL;
        }

        inline bool close() {
            bool success = CloseHandle(hEvent);
			hEvent = NULL;
			return success;
        }

        inline bool wait(DWORD timeout_millis = Infinite) {
            if (!hEvent) return true;
            return WaitForSingleObject(hEvent, timeout_millis) != WAIT_TIMEOUT;
        }

        inline bool signal() {
            return SetEvent(hEvent) != 0;
        }

        inline bool reset() {
            return ResetEvent(hEvent) != 0;
        }

        inline operator bool() const {
            return hEvent != NULL;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global event functions
	
	inline void reset(Event &e) {
		e.reset();
	}

    inline bool waitForAll(Event eventList[], DWORD count, DWORD timeout_millis)
	{
		return WaitForMultipleObjects(count, &eventList[0].hEvent, true, timeout_millis) != WAIT_TIMEOUT;
	}
    
    inline bool waitForAll(Event eventList[], DWORD count)
	{
		return waitForAll(eventList, count, Event::Infinite);
	}
    
	inline int waitForAny(Event eventList[], DWORD count, DWORD timeout_millis)
	{
        DWORD result = WaitForMultipleObjects(count, &eventList[0].hEvent, false, timeout_millis);

        if (/*result >= WAIT_OBJECT_0 && */result < WAIT_OBJECT_0 + count)
    		return result - WAIT_OBJECT_0;

        return -1;
	}
    
	inline int waitForAny(Event eventList[], DWORD count)
	{
		return waitForAny(eventList, count, Event::Infinite);
	}
    
#else

    struct Event: public Mutex
    {
        typedef pthread_cond_t* Handle;
		enum { Infinite = LONG_MAX };
        pthread_cond_t data, *hEvent;

        Event():
            hEvent(NULL) {}

        Event(bool initial) {
            SEQAN_DO_SYS(open(initial));
        }

        ~Event() {
            if (*this)
                SEQAN_DO_SYS(close());
        }

		Event(Event const &origin):
			Mutex()
		{
            // resource sharing is not yet supported (performance reason)
            // it needs a reference counting technique
            if (origin) {
                data = origin.data;
                const_cast<Event&>(origin).hEvent = NULL;
                hEvent = &data;
            } else
                hEvent = NULL;
        }

        inline Event& operator=(Event const &origin) {
            // resource sharing is not yet supported (performance reason)
            // it needs a reference counting technique
            if (origin) {
                data = origin.data;
                const_cast<Event&>(origin).hEvent = NULL;
                hEvent = &data;
            } else
                hEvent = NULL;
            return *this;
        }

        inline bool open(bool initial = false)
        {
            if (Mutex::open() && !pthread_cond_init(&data, NULL) && (hEvent = &data)) {
                if (initial) return signal();
                return true;
            } else
                return false;
        }

        inline bool close() {
            bool success = (pthread_cond_destroy(hEvent) == 0);
		    success &= Mutex::close();
            hEvent = NULL;
            return success;
		}

        inline bool wait() {
            if (!hEvent) return true;
            Mutex::lock();
            return !pthread_cond_wait(hEvent, Mutex::hMutex);
        }

        inline bool wait(long timeout_millis) {
            if (timeout_millis != Infinite) {
                timespec ts;
                ts.tv_sec = timeout_millis / 1000;
                ts.tv_nsec = (timeout_millis % 1000) * 1000;
                Mutex::lock();
				return pthread_cond_timedwait(hEvent, Mutex::hMutex, &ts) != ETIMEDOUT;
            } else
                return wait();
        }

        inline bool signal() {
            return !pthread_cond_broadcast(hEvent);
        }

        inline operator bool() const {
            return hEvent != NULL;
        }
    };
    
#endif


	//////////////////////////////////////////////////////////////////////////////
	// global event functions

	inline bool open(Event &e, bool initial) {
		return e.open(initial);
	}

	inline bool open(Event &e) {
		return open(e, false);
	}

	inline bool close(Event &e) {
		return e.close();
	}

	inline bool waitFor(Event &e) {
		return e.wait();
	}

    template < typename TTime >
	inline bool waitFor(Event &e, TTime timeout_millis) {
        #ifdef disabledSEQAN_PROFILE
			double begin = sysTime();
			bool b = e.wait(timeout_millis);
			double end = sysTime();
            if (begin != end)
                ::std::cerr << "waitTime: " << end - begin << ::std::endl;
			return b;
        #else
            return e.wait(timeout_millis);
        #endif
	}

	inline bool signal(Event &e) {
		return e.signal();
	}

/*
	//////////////////////////////////////////////////////////////////////////////
	// emulate events in a singlethreaded environment

	struct DummyEvent {
		typedef	void Handle;
		DummyEvent(bool initial = false) {}
		inline bool wait(unsigned timeOut = NULL) { return true; }
		inline void reset() {}
		inline void signal() {}
	};

	//////////////////////////////////////////////////////////////////////////////
	// global dummy event functions

	template < typename TCount >
	inline bool waitForAll(DummyEvent eventList[], TCount count) {
		return true;
	}

	template < typename TCount, typename TTime >
	inline bool waitForAll(DummyEvent eventList[], TCount count, TTime timeout_millis) {
		return true;
	}

	template < typename TCount >
	inline TCount waitForAny(DummyEvent eventList[], TCount count) {
		return 0;
	}
	template < typename TCount, typename TTime >
	inline TCount waitForAny(DummyEvent eventList[], TCount count, TTime timeout_millis) {
		return 0;
	}
*/
}

#endif
