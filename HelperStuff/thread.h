#ifndef _thread_h
#define _thread_h

/*=========================================================================
                                                                                
  Program:   thread
  Module:    $RCSfile: thread.h,v $
  Language:  C++
  Date:      $Date: 2014/04/06 14:14:50 $
  Version:   $Revision: 1.0 $
                                                                                
  Copyright (c) BME
   
                                                                                
=========================================================================*/




#include <process.h>    /* _beginthread, _endthread */
#include "Timer.h"
//#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <iostream>
#include <string>
#include <cstdarg>
//#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#pragma once

/*!
 ... text ...
*/


 /// \todo More Documentation

/*!
 * \brief Timer + CriticalSection objects
 *        Base Timer is used to synchronize all threads !
 * the critical objects are used as wait signals at times thad shared opjects are filled 
 *
 * Remarks
 * The threads of a single process can use a critical section object for mutual-exclusion synchronization. 
 * There is no guarantee about the order in which threads will obtain ownership of the critical section, 
 * however, the system will be fair to all threads.
 * @param   
 */


class SyncTimerObject:public Timer
{
public:
	SyncTimerObject()
	{
		InitializeCriticalSection(&m_CriticalSection1);
		InitializeCriticalSection(&m_CriticalSection2);
		InitializeCriticalSection(&m_CriticalSection3);
	}

	virtual ~SyncTimerObject()
	{
		// Release resources used by the critical section object.
		DeleteCriticalSection(&m_CriticalSection1);
		DeleteCriticalSection(&m_CriticalSection2);
		DeleteCriticalSection(&m_CriticalSection3);
	}

/*!
 * \brief StartObjectWriteSyncronisationSection
 *  
 * EnterCriticalSection:: Waits for ownership of the specified critical section object. The function returns when the calling thread is granted ownership.
 *
 * Remarks
 * The threads of a single process can use a critical section object for mutual-exclusion synchronization. 
 * The process is responsible for allocating the memory used by a critical section object, 
 * @param   
 */
	void StartObjectWriteSyncronisationSection(int ObjectNR) 
	{

		
		switch(ObjectNR)
		{
		case 1:
			EnterCriticalSection( &m_CriticalSection1 );  
			break;
		case 2:
			EnterCriticalSection( &m_CriticalSection2 );  
			break;
		case 3:
			EnterCriticalSection( &m_CriticalSection3 );  
			break;
		}
		
	}

/*!
 * \brief StopObjectWriteSyncronisationSection
 *  
 * LeaveCriticalSection:: Releases ownership of the specified critical section object.
 *
 * Remarks
 * The threads of a single process can use a critical-section object for mutual-exclusion synchronization. 
 * The process is responsible for allocating the memory used by a critical-section object
 * @param   
 */
	
	void StopObjectWriteSyncronisationSection(int ObjectNR) 
	{
		switch(ObjectNR)
		{
		case 1:
			LeaveCriticalSection( &m_CriticalSection1 );  
			break;
		case 2:
			LeaveCriticalSection( &m_CriticalSection2 );  
			break;
		case 3:
			LeaveCriticalSection( &m_CriticalSection3 );  
			break;
		}
	}

/*!
 * \brief WaitForObjectWriteSyncronisationSection
 *  
 *  Waits for ownership of the specified critical section object and directly Releases ownership of the specified critical section object .
 *
 * Remarks
 * calling both EnterCriticalSection and LeaveCriticalSection
 * 
 
 */
	 void WaitForObjectWriteSyncronisationSection(int ObjectNR)
	 {
		 StartObjectWriteSyncronisationSection( ObjectNR);
		 StopObjectWriteSyncronisationSection( ObjectNR);
	 }

private:
	CRITICAL_SECTION m_CriticalSection1; 
	CRITICAL_SECTION m_CriticalSection2; 
	CRITICAL_SECTION m_CriticalSection3; 
public:

};


/*!
 * \brief thread Base Class for Multithreading objects
 *    This is a templated class with Input and Output  
 * The timing/synchronisation object is passed by the function SetSync
 *
 * Remarks
 *
 *
 * @param   
 */


template <class Tin,class Tout> class thread{
private:
	HANDLE threadhandle;
	unsigned threadid;

public:
/*!
 ... text ...
*/
	thread();
	~thread();
	
	static thread<Tin,Tout> * New( thread<Tin,Tout> * );
	virtual void SetSync( SyncTimerObject * SyncOpject );
/*!
 ... SetInput ...
*/
	virtual void SetInput( Tin * InputData );
	virtual void SetOutput( Tout * OutputData );
	virtual void Update(); //Prepare thread
	Tout *  GetOutput();
	
	Tin *m_InputData;
	Tout *m_OutputData;
	SyncTimerObject *m_SyncObject;

public:
	virtual void Initialize(){};
/*
va_list

Type to hold information about variable arguments

This type is used as a parameter for the macros defined in <cstdarg> to retrieve the additional arguments of a function.

va_start initializes an object of this type in such a way that subsequent calls to va_arg sequentially retrieve
         the additional arguments passed to the function.

Before a function that has initialized a va_list object with va_start returns, 
   the va_end macro shall be invoked.
   */
/*!
void Initialize(const char *fmt, ...)
Pass a variable parameter list 
the format string gives a format code for each parameter  <BR>
 - W = word  <BR>
 - D = double  <BR>
*/

	virtual void Initialize(const char *fmt, ...)
	{
		va_list args;
		va_start(args, fmt);

		while (*fmt != '\0') {
			if (*fmt == 'd') {
				int i = va_arg(args, int);
				std::cout << i << '\n';
			} else if (*fmt == 's') {
				char * s = va_arg(args, char*);
				std::cout << s << '\n';
			}
			++fmt;
		}
		va_end(args);
	};
	//void BeginThread(thread<Tin,Tout> *_Thread)

/*!
    Creates a thread that begins execution on _BeginThread
*/
	void _BeginThread()

	{
		threadhandle = (HANDLE)_beginthreadex(NULL, 0,thread<Tin,Tout>::ThreadStaticEntryPoint,this, CREATE_SUSPENDED, &threadid);
	}
/*!
    Start/Resume  the execution of the thread .
*/
	void _ResumeThread(void){ResumeThread( threadhandle );};
/*!
    Waits until the specified thread is ready
*/	
	void _WaitForSingleObject(void){WaitForSingleObject(threadhandle, INFINITE);};
	
	static unsigned __stdcall ThreadStaticEntryPoint(void * pThis)
	{
		thread * pthX = (thread*)pThis;		// the tricky cast

		pthX->ThreadEntryPoint();			// now call the true entry-point-function

		// A thread terminates automatically if it completes execution,
		// or it can terminate itself with a call to _endthread().
		return 1;          // the thread exit code
	}

	virtual void ThreadEntryPoint()
	
	{

	};
		// This is the desired entry-point-function but to get
		// here we have to use a 2 step procedure involving
		// the ThreadStaticEntryPoint() function.
	

	static thread *m_thread;
 	private:
 
	//Timer *m_timer;
	const char *m_name;
};
template <class Tin,class Tout>  thread<Tin,Tout>::thread(void)
{
	m_InputData=NULL;
	m_OutputData=NULL;
};
template  <class Tin,class Tout>  thread<Tin,Tout>::~thread(void)
{

}

/*!
	The Base class is Templated so any type of class can be used as input <BR>
	This Class is used to Pass the input data to the thread

	the class can be used as follows
	m_InputData->SharedObjects::SetTimeVideo(&_timeVideo);

*/
template  <class Tin,class Tout> void thread<Tin,Tout>::SetInput(  Tin  * InputData )
{
	m_InputData=InputData;
}



/*!
For synchronization we use the SyncTimerObject  <BR>
this object exist of a Timer and a number of CRITICAL_SECTION objects  <BR>
with those CRITICAL_SECTION objects  we can let on thread wait until another thread is ready filling a data block <BR>

*/
template  <class Tin,class Tout> void thread<Tin,Tout>::SetSync(  SyncTimerObject *SharedObjects)
{
	m_SyncObject=SharedObjects;
}

template  <class Tin,class Tout> Tout  *  thread<Tin,Tout>::GetOutput( )
{
	return m_OutputData;
}


/*!
	The Base class is Templated so any type of class can be used as Output, <BR>
	
	the class can be used as follows
	
	m_OutputData->TiepieObject::SetTimeCh1(_timeCh1);

*/

template  <class Tin,class Tout> void thread<Tin,Tout>::SetOutput(  Tout  * OutputData )
{
	m_OutputData=OutputData;
}
template <class Tin,class Tout> void thread<Tin,Tout>::Update(void)
{

};

//template <class Tin,class Tout>  thread<Tin,Tout> *  thread<Tin,Tout>::m_thread=0;

template  <class Tin,class Tout>  thread<Tin,Tout> *  thread<Tin,Tout>::New( thread<Tin,Tout> * _thread=NULL)
{
	if(m_thread==NULL)
	{
		m_thread= new thread;
	}
	return m_thread;
}

#endif