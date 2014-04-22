#ifndef __TIEPIE_H__
#define __TIEPIE_H__

/*
 * This header file contains the declarations for the low level routines
 * for the TiePie engineering instrument dll's
 */

#include <windows.h>


/*
 * declarations to match pascal-types to c-types
 */
typedef unsigned char      byte;
typedef unsigned short int word;
typedef unsigned int       dword;

typedef signed   long int  integer;
typedef signed long long   int64;

typedef void *             pointer;

/*
 * The instruments
 */
typedef enum {
    dtHandyprobe2,
    dtHandyscope2,
    dtHandyscope3,
    dtHandyscope4,
    dtHS508,
    dtHS801,
    dtTP112,
    dtTP208,
    dtTP508,
    dtTP801,
    dtPCI801,
    dtMultiCh
} TDeviceTypes;

// Return codes = bit values. Non-zero is an error:
#define E_NO_ERRORS              0x0000
#define E_NO_HARDWARE            0x0001
#define E_NOT_INITIALIZED        0x0002
#define E_NOT_SUPPORTED          0x0004
#define E_NO_GENERATOR           0x0008
#define E_INVALID_CHANNEL        0x0010
#define E_INVALID_VALUE          0x0020
#define E_I2C_ERROR              0x0040
#define E_I2C_INVALID_ADDRESS    0x0080
#define E_I2C_INVALID_SIZE       0x0100
#define E_I2C_NO_ACKNOWLEDGE     0x0200

// Channel constants:
#define lCh1                      1
#define lCh2                      2
#define lCh3                      3
#define lCh4                      4

// Measure Mode constants, do NOT use for anything else than G/SetMeasureMode().
// When measuring more channels, simply add the required constants together : e.g. mmCh1 + mmCh2 to measure both Ch1 and Ch2
#define mmCh1                     1
#define mmCh2                     2
#define mmCh3                     4
#define mmCh4                     8
#define mmCh12                    3

// Coupling constants:
#define lctAC                     0
#define lctDC                     1

// Trigger mode constants:
#define ltkRising                 0 // Trigger when signal passes level in positive direction
#define ltkFalling                1 // Trigger when signal passes level in negative direction
#define ltkInWindow               2 // Trigger when signal enters marked window
#define ltkOutWindow              3 // Trigger when signal leaves marked window
#define ltkTVLine                 4 // Trigger at TV line sync pulse
#define ltkTVFieldOdd             5 // Trigger at Odd TV field sync pulse
#define ltkTVFieldEven            6 // Trigger at Even TV field sync pulse
#define ltkNet                    7 // Trigger on peaks/dips in the line net
#define ltkNthPulse               8 // Trigger on the nth pulse after start of measuring
#define ltkPulseLargerM           9 // Trigger if a pulse is larger than M
#define ltkPulseSmallerM         10 // Trigger if a pulse is smaller than M
#define ltkGlitch                11 // Trigger when the change between two samples is larger than x

// Trigger source constants:
#define ltsCh1                    0 // Trigger on Ch1
#define ltsCh2                    1 // Trigger on Ch2
#define ltsCh3                    2 // Trigger on Ch3
#define ltsCh4                    3 // Trigger on Ch4
#define ltsExternal               4 // Trigger on digital external trigger input
#define ltsAnalogExt              5 // Trigger on analog external trigger input
#define ltsAnd                    6 // Trigger on Ch1 AND Ch2
#define ltsOr                     7 // Trigger on Ch1 OR  Ch2
#define ltsXor                    8 // Trigger on Ch1 XOR Ch2
#define ltsNoTrig                 9 // Do not wait for trigger, measure immediately
#define ltsNeverTrig             10 // not used
#define ltsPxiExt                11 // Trigger on external PXI trigger inputs
#define ltsGenStart              12 // Trigger on Generator Start Signal
#define ltsGenStop               13 // Trigger on Generator Stop Signal
#define ltsGenNew                14 // Trigger on Generator New Period
#define ltsUnknown               15 // not used

// Transfer mode constants:
#define ltmBlock                  0
#define ltmStream                 1

// Instrument Configuration constants:
#define licHS3Norm                0
#define licHS3256K                1
#define licHS3512K                2

// FuncGen Modes:
#define lfmDDS                    1
#define lfmLinear                 2

// Signal types
#define lstSine                   0
#define lstTriangle               1
#define lstSquare                 2
#define lstDC                     3
#define lstNoise                  4
#define lstArbitrary              5


/*
 * definition of Callback function
 */
typedef void (*TDataReady) (void);


/*
 * Properly decorate the API functions for either import (default) or export.
 */
#ifdef BUILD_TIEPIE
#   define TP_STORAGE_CLASS DLLEXPORT
#else
#   define TP_STORAGE_CLASS DLLIMPORT
#endif

/*
 *  Adjust naming convention to match language in use.
 */
#ifdef __cplusplus
#   define TP_EXTERNC extern "C"
#else
#   define TP_EXTERNC extern
#endif

/*
 * Use the proper storage-class attributes (and order) for the compiler
 * in use.
 */
#if (defined(_MSC_VER) || (__BORLANDC__ >= 0x0550) || \
	defined(__LCC__) || defined(__WATCOMC__) || \
	(defined(__GNUC__) && defined(__declspec)))
#   define DLLIMPORT __declspec(dllimport)
#   define DLLEXPORT __declspec(dllexport)
#   define TP_EXTERN(RTYPE) TP_EXTERNC TP_STORAGE_CLASS RTYPE WINAPI
#elif defined(__BORLANDC__)
    /* Pre-5.5 Borland requires the attributes be placed after the */
    /* return type instead. */
#   define DLLIMPORT __import
#   define DLLEXPORT __export
#   define TP_EXTERN(RTYPE) TP_EXTERNC RTYPE TP_STORAGE_CLASS WINAPI
#else
#   define DLLIMPORT
#   define DLLEXPORT
#   define TP_EXTERN(RTYPE) TP_EXTERNC RTYPE WINAPI
#endif


/*
 * ============================================================================
 * Declarations for functions from dynamically loaded librarie(s)
 * ============================================================================
 */

// Open / Close the instrument:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_InitInstrument) ( word wAddress );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_InitInstrument InitInstrument;
  #endif
#else
    TP_EXTERN( word ) InitInstrument( word wAddress );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_ExitInstrument) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_ExitInstrument ExitInstrument;
  #endif
#else
    TP_EXTERN( word ) ExitInstrument( void );
#endif
 
 
// Get information about the instrument:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetCalibrationDate) ( dword* dwDate );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetCalibrationDate GetCalibrationDate;
  #endif
#else
    TP_EXTERN( word ) GetCalibrationDate( dword* dwDate );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetSerialNumber) ( dword* dwSerialNumber );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetSerialNumber GetSerialNumber;
  #endif
#else
    TP_EXTERN( word ) GetSerialNumber( dword* dwSerialNumber );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetAvailableSensitivities) ( double* dSensitivities );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetAvailableSensitivities GetAvailableSensitivities;
  #endif
#else
    TP_EXTERN( word ) GetAvailableSensitivities( double* dSensitivities );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetAvailableResolutions) ( double* dResolutions );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetAvailableResolutions GetAvailableResolutions;
  #endif
#else
    TP_EXTERN( word ) GetAvailableResolutions( double* dResolutions );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetNrChannels) ( word* wValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetNrChannels GetNrChannels;
  #endif
#else
    TP_EXTERN( word ) GetNrChannels( word* wValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef double (WINAPI * LPFN_GetMaxSampleFrequencyF) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetMaxSampleFrequencyF GetMaxSampleFrequencyF;
  #endif
#else
    TP_EXTERN( double ) GetMaxSampleFrequencyF( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef dword (WINAPI * LPFN_GetMaxRecordLength) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetMaxRecordLength GetMaxRecordLength;
  #endif
#else
    TP_EXTERN( dword ) GetMaxRecordLength( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef dword (WINAPI * LPFN_GetDCLevelStatus) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetDCLevelStatus GetDCLevelStatus;
  #endif
#else
    TP_EXTERN( dword ) GetDCLevelStatus( void );
#endif
 
 
// Measurement routines:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_ADC_Start) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_ADC_Start ADC_Start;
  #endif
#else
    TP_EXTERN( word ) ADC_Start( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_ADC_Running) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_ADC_Running ADC_Running;
  #endif
#else
    TP_EXTERN( word ) ADC_Running( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_ADC_Abort) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_ADC_Abort ADC_Abort;
  #endif
#else
    TP_EXTERN( word ) ADC_Abort( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_ADC_ForceTrig) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_ADC_ForceTrig ADC_ForceTrig;
  #endif
#else
    TP_EXTERN( word ) ADC_ForceTrig( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_ADC_Triggered) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_ADC_Triggered ADC_Triggered;
  #endif
#else
    TP_EXTERN( word ) ADC_Triggered( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_ADC_Ready) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_ADC_Ready ADC_Ready;
  #endif
#else
    TP_EXTERN( word ) ADC_Ready( void );
#endif
 
 
// Retrieve the measured data:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_ADC_GetDataCh) ( word wCh , word* wChData );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_ADC_GetDataCh ADC_GetDataCh;
  #endif
#else
    TP_EXTERN( word ) ADC_GetDataCh( word wCh , word* wChData );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_ADC_GetDataVoltCh) ( word wCh , double* dChData );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_ADC_GetDataVoltCh ADC_GetDataVoltCh;
  #endif
#else
    TP_EXTERN( word ) ADC_GetDataVoltCh( word wCh , double* dChData );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetDigitalInputValues) ( word* wValues );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetDigitalInputValues GetDigitalInputValues;
  #endif
#else
    TP_EXTERN( word ) GetDigitalInputValues( word* wValues );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetOneDigitalValue) ( dword dwIndex , word* wValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetOneDigitalValue GetOneDigitalValue;
  #endif
#else
    TP_EXTERN( word ) GetOneDigitalValue( dword dwIndex , word* wValue );
#endif
 
 
// Streaming measurements:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetDataReadyCallback) ( TDataReady pAddress );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetDataReadyCallback SetDataReadyCallback;
  #endif
#else
    TP_EXTERN( word ) SetDataReadyCallback( TDataReady pAddress );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetDataReadyEvent) ( HANDLE hEvent );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetDataReadyEvent SetDataReadyEvent;
  #endif
#else
    TP_EXTERN( word ) SetDataReadyEvent( HANDLE hEvent );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetTransferMode) ( dword dwMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetTransferMode SetTransferMode;
  #endif
#else
    TP_EXTERN( word ) SetTransferMode( dword dwMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetTransferMode) ( dword* dwMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetTransferMode GetTransferMode;
  #endif
#else
    TP_EXTERN( word ) GetTransferMode( dword* dwMode );
#endif
 
 
// Control the input resolution in bits:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetResolution) ( byte byValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetResolution SetResolution;
  #endif
#else
    TP_EXTERN( word ) SetResolution( byte byValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetResolution) ( byte* byValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetResolution GetResolution;
  #endif
#else
    TP_EXTERN( word ) GetResolution( byte* byValue );
#endif
 
 
// Control the instrument configuration:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetInstrumentConfig) ( word wMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetInstrumentConfig SetInstrumentConfig;
  #endif
#else
    TP_EXTERN( word ) SetInstrumentConfig( word wMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetInstrumentConfig) ( word* wMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetInstrumentConfig GetInstrumentConfig;
  #endif
#else
    TP_EXTERN( word ) GetInstrumentConfig( word* wMode );
#endif
 
 
// Control which channels are measured:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetMeasureMode) ( byte* byMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetMeasureMode GetMeasureMode;
  #endif
#else
    TP_EXTERN( word ) GetMeasureMode( byte* byMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetMeasureMode) ( byte byMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetMeasureMode SetMeasureMode;
  #endif
#else
    TP_EXTERN( word ) SetMeasureMode( byte byMode );
#endif
 
 
// Control the time base:
#ifdef TP_TYPEDEFS
    typedef dword (WINAPI * LPFN_GetRecordLength) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetRecordLength GetRecordLength;
  #endif
#else
    TP_EXTERN( dword ) GetRecordLength( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetRecordLength) ( dword dwRecordLength );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetRecordLength SetRecordLength;
  #endif
#else
    TP_EXTERN( word ) SetRecordLength( dword dwRecordLength );
#endif
 
#ifdef TP_TYPEDEFS
    typedef dword (WINAPI * LPFN_GetPostSamples) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetPostSamples GetPostSamples;
  #endif
#else
    TP_EXTERN( dword ) GetPostSamples( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetPostSamples) ( dword dwPostSampleCount );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetPostSamples SetPostSamples;
  #endif
#else
    TP_EXTERN( word ) SetPostSamples( dword dwPostSampleCount );
#endif
 
#ifdef TP_TYPEDEFS
    typedef double (WINAPI * LPFN_GetSampleFrequencyF) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetSampleFrequencyF GetSampleFrequencyF;
  #endif
#else
    TP_EXTERN( double ) GetSampleFrequencyF( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetSampleFrequencyF) ( double* dFreq );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetSampleFrequencyF SetSampleFrequencyF;
  #endif
#else
    TP_EXTERN( word ) SetSampleFrequencyF( double* dFreq );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetExternalClock) ( word* wMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetExternalClock GetExternalClock;
  #endif
#else
    TP_EXTERN( word ) GetExternalClock( word* wMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetExternalClock) ( word wMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetExternalClock SetExternalClock;
  #endif
#else
    TP_EXTERN( word ) SetExternalClock( word wMode );
#endif
 
 
// Control the analog input channels:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetSensitivity) ( byte byCh , double* dSensitivity );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetSensitivity GetSensitivity;
  #endif
#else
    TP_EXTERN( word ) GetSensitivity( byte byCh , double* dSensitivity );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetSensitivity) ( byte byCh , double* dSensitivity );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetSensitivity SetSensitivity;
  #endif
#else
    TP_EXTERN( word ) SetSensitivity( byte byCh , double* dSensitivity );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetAutoRanging) ( byte byCh , byte* byMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetAutoRanging GetAutoRanging;
  #endif
#else
    TP_EXTERN( word ) GetAutoRanging( byte byCh , byte* byMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetAutoRanging) ( byte byCh , byte byMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetAutoRanging SetAutoRanging;
  #endif
#else
    TP_EXTERN( word ) SetAutoRanging( byte byCh , byte byMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetCoupling) ( byte byCh , byte* byMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetCoupling GetCoupling;
  #endif
#else
    TP_EXTERN( word ) GetCoupling( byte byCh , byte* byMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetCoupling) ( byte byCh , byte byMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetCoupling SetCoupling;
  #endif
#else
    TP_EXTERN( word ) SetCoupling( byte byCh , byte byMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetDcLevel) ( byte byCh , double* dDcLevel );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetDcLevel GetDcLevel;
  #endif
#else
    TP_EXTERN( word ) GetDcLevel( byte byCh , double* dDcLevel );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetDcLevel) ( byte byCh , double dDcLevel );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetDcLevel SetDcLevel;
  #endif
#else
    TP_EXTERN( word ) SetDcLevel( byte byCh , double dDcLevel );
#endif
 
 
// Control the trigger system:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetTriggerSource) ( byte* byTrigSource );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetTriggerSource GetTriggerSource;
  #endif
#else
    TP_EXTERN( word ) GetTriggerSource( byte* byTrigSource );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetTriggerSource) ( byte byTrigSource );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetTriggerSource SetTriggerSource;
  #endif
#else
    TP_EXTERN( word ) SetTriggerSource( byte byTrigSource );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetTriggerMode) ( byte* byMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetTriggerMode GetTriggerMode;
  #endif
#else
    TP_EXTERN( word ) GetTriggerMode( byte* byMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetTriggerMode) ( byte byMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetTriggerMode SetTriggerMode;
  #endif
#else
    TP_EXTERN( word ) SetTriggerMode( byte byMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetTriggerModeCh) ( byte byCh , byte* byMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetTriggerModeCh GetTriggerModeCh;
  #endif
#else
    TP_EXTERN( word ) GetTriggerModeCh( byte byCh , byte* byMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetTriggerModeCh) ( byte byCh , byte byMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetTriggerModeCh SetTriggerModeCh;
  #endif
#else
    TP_EXTERN( word ) SetTriggerModeCh( byte byCh , byte byMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetTriggerLevel) ( byte byCh , double* dLevel );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetTriggerLevel GetTriggerLevel;
  #endif
#else
    TP_EXTERN( word ) GetTriggerLevel( byte byCh , double* dLevel );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetTriggerLevel) ( byte byCh , double dLevel );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetTriggerLevel SetTriggerLevel;
  #endif
#else
    TP_EXTERN( word ) SetTriggerLevel( byte byCh , double dLevel );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetTriggerHys) ( byte byCh , double* dHysteresis );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetTriggerHys GetTriggerHys;
  #endif
#else
    TP_EXTERN( word ) GetTriggerHys( byte byCh , double* dHysteresis );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetTriggerHys) ( byte byCh , double dHysteresis );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetTriggerHys SetTriggerHys;
  #endif
#else
    TP_EXTERN( word ) SetTriggerHys( byte byCh , double dHysteresis );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetPXITriggerEnables) ( byte* byEnables );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetPXITriggerEnables GetPXITriggerEnables;
  #endif
#else
    TP_EXTERN( word ) GetPXITriggerEnables( byte* byEnables );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetPXITriggerEnables) ( byte byEnables );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetPXITriggerEnables SetPXITriggerEnables;
  #endif
#else
    TP_EXTERN( word ) SetPXITriggerEnables( byte byEnables );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetPXITriggerSlopes) ( byte* bySlopes );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetPXITriggerSlopes GetPXITriggerSlopes;
  #endif
#else
    TP_EXTERN( word ) GetPXITriggerSlopes( byte* bySlopes );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetPXITriggerSlopes) ( byte bySlopes );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetPXITriggerSlopes SetPXITriggerSlopes;
  #endif
#else
    TP_EXTERN( word ) SetPXITriggerSlopes( byte bySlopes );
#endif
 
 
// Control the digital outputs:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetDigitalOutputs) ( byte byValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetDigitalOutputs SetDigitalOutputs;
  #endif
#else
    TP_EXTERN( word ) SetDigitalOutputs( byte byValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetDigitalOutputs) ( byte* byValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetDigitalOutputs GetDigitalOutputs;
  #endif
#else
    TP_EXTERN( word ) GetDigitalOutputs( byte* byValue );
#endif
 
 
// Square Wave Generator:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetSquareWaveGenStatus) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetSquareWaveGenStatus GetSquareWaveGenStatus;
  #endif
#else
    TP_EXTERN( word ) GetSquareWaveGenStatus( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef dword (WINAPI * LPFN_GetSquareWaveGenFrequency) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetSquareWaveGenFrequency GetSquareWaveGenFrequency;
  #endif
#else
    TP_EXTERN( dword ) GetSquareWaveGenFrequency( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetSquareWaveGenFrequency) ( dword* dwFreq );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetSquareWaveGenFrequency SetSquareWaveGenFrequency;
  #endif
#else
    TP_EXTERN( word ) SetSquareWaveGenFrequency( dword* dwFreq );
#endif
 
 
// Get information about the function generator:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFunctionGenStatus) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFunctionGenStatus GetFunctionGenStatus;
  #endif
#else
    TP_EXTERN( word ) GetFunctionGenStatus( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFuncGenMaxAmplitude) ( double* dMaxAmplitude );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFuncGenMaxAmplitude GetFuncGenMaxAmplitude;
  #endif
#else
    TP_EXTERN( word ) GetFuncGenMaxAmplitude( double* dMaxAmplitude );
#endif
 
 
// Control the Arbitrary Waveform Generator:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFuncGenOutputOn) ( word* wValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFuncGenOutputOn GetFuncGenOutputOn;
  #endif
#else
    TP_EXTERN( word ) GetFuncGenOutputOn( word* wValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetFuncGenOutputOn) ( word wValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetFuncGenOutputOn SetFuncGenOutputOn;
  #endif
#else
    TP_EXTERN( word ) SetFuncGenOutputOn( word wValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFuncGenEnable) ( word* wValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFuncGenEnable GetFuncGenEnable;
  #endif
#else
    TP_EXTERN( word ) GetFuncGenEnable( word* wValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetFuncGenEnable) ( word wValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetFuncGenEnable SetFuncGenEnable;
  #endif
#else
    TP_EXTERN( word ) SetFuncGenEnable( word wValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFuncGenSignalType) ( word* wValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFuncGenSignalType GetFuncGenSignalType;
  #endif
#else
    TP_EXTERN( word ) GetFuncGenSignalType( word* wValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetFuncGenSignalType) ( word wValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetFuncGenSignalType SetFuncGenSignalType;
  #endif
#else
    TP_EXTERN( word ) SetFuncGenSignalType( word wValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetFuncGenMode) ( dword dwMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetFuncGenMode SetFuncGenMode;
  #endif
#else
    TP_EXTERN( word ) SetFuncGenMode( dword dwMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFuncGenMode) ( dword* dwMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFuncGenMode GetFuncGenMode;
  #endif
#else
    TP_EXTERN( word ) GetFuncGenMode( dword* dwMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFuncGenAmplitude) ( double* dValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFuncGenAmplitude GetFuncGenAmplitude;
  #endif
#else
    TP_EXTERN( word ) GetFuncGenAmplitude( double* dValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetFuncGenAmplitude) ( double dValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetFuncGenAmplitude SetFuncGenAmplitude;
  #endif
#else
    TP_EXTERN( word ) SetFuncGenAmplitude( double dValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFuncGenAmplitudeRange) ( byte* byValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFuncGenAmplitudeRange GetFuncGenAmplitudeRange;
  #endif
#else
    TP_EXTERN( word ) GetFuncGenAmplitudeRange( byte* byValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetFuncGenAmplitudeRange) ( byte byValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetFuncGenAmplitudeRange SetFuncGenAmplitudeRange;
  #endif
#else
    TP_EXTERN( word ) SetFuncGenAmplitudeRange( byte byValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFuncGenDCOffset) ( double* dValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFuncGenDCOffset GetFuncGenDCOffset;
  #endif
#else
    TP_EXTERN( word ) GetFuncGenDCOffset( double* dValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetFuncGenDCOffset) ( double dValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetFuncGenDCOffset SetFuncGenDCOffset;
  #endif
#else
    TP_EXTERN( word ) SetFuncGenDCOffset( double dValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFuncGenSymmetry) ( double* dValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFuncGenSymmetry GetFuncGenSymmetry;
  #endif
#else
    TP_EXTERN( word ) GetFuncGenSymmetry( double* dValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetFuncGenSymmetry) ( double dValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetFuncGenSymmetry SetFuncGenSymmetry;
  #endif
#else
    TP_EXTERN( word ) SetFuncGenSymmetry( double dValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFuncGenFrequency) ( double* dValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFuncGenFrequency GetFuncGenFrequency;
  #endif
#else
    TP_EXTERN( word ) GetFuncGenFrequency( double* dValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetFuncGenFrequency) ( double* dValue );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetFuncGenFrequency SetFuncGenFrequency;
  #endif
#else
    TP_EXTERN( word ) SetFuncGenFrequency( double* dValue );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetFuncGenTrigSource) ( byte byTrigSource );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetFuncGenTrigSource SetFuncGenTrigSource;
  #endif
#else
    TP_EXTERN( word ) SetFuncGenTrigSource( byte byTrigSource );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetFuncGenTrigSource) ( byte* byTrigSource );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetFuncGenTrigSource GetFuncGenTrigSource;
  #endif
#else
    TP_EXTERN( word ) GetFuncGenTrigSource( byte* byTrigSource );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_FillFuncGenMemory) ( dword dwSampleCount , word* FuncGenData );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_FillFuncGenMemory FillFuncGenMemory;
  #endif
#else
    TP_EXTERN( word ) FillFuncGenMemory( dword dwSampleCount , word* FuncGenData );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_FuncGenBurst) ( word wNrPeriods );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_FuncGenBurst FuncGenBurst;
  #endif
#else
    TP_EXTERN( word ) FuncGenBurst( word wNrPeriods );
#endif
 
 
// Ohm measurements ( not available by default ):
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetupOhmMeasurements) ( word wMode );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetupOhmMeasurements SetupOhmMeasurements;
  #endif
#else
    TP_EXTERN( word ) SetupOhmMeasurements( word wMode );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetOhmValues) ( double* dValue1 , double* dValue2 );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetOhmValues GetOhmValues;
  #endif
#else
    TP_EXTERN( word ) GetOhmValues( double* dValue1 , double* dValue2 );
#endif
 
 
// I2C routines:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_I2CWrite) ( dword dwAddress , pointer pBuf , dword dwSize );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_I2CWrite I2CWrite;
  #endif
#else
    TP_EXTERN( word ) I2CWrite( dword dwAddress , pointer pBuf , dword dwSize );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_I2CWriteNoStop) ( dword dwAddress , pointer pBuf , dword dwSize );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_I2CWriteNoStop I2CWriteNoStop;
  #endif
#else
    TP_EXTERN( word ) I2CWriteNoStop( dword dwAddress , pointer pBuf , dword dwSize );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_I2CRead) ( dword dwAddress , pointer pBuf , dword dwSize );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_I2CRead I2CRead;
  #endif
#else
    TP_EXTERN( word ) I2CRead( dword dwAddress , pointer pBuf , dword dwSize );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_I2CReadNoStop) ( dword dwAddress , pointer pBuf , dword dwSize );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_I2CReadNoStop I2CReadNoStop;
  #endif
#else
    TP_EXTERN( word ) I2CReadNoStop( dword dwAddress , pointer pBuf , dword dwSize );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_I2CGetSpeed) ( dword* dwSpeed );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_I2CGetSpeed I2CGetSpeed;
  #endif
#else
    TP_EXTERN( word ) I2CGetSpeed( dword* dwSpeed );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_I2CSetSpeed) ( dword* dwSpeed );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_I2CSetSpeed I2CSetSpeed;
  #endif
#else
    TP_EXTERN( word ) I2CSetSpeed( dword* dwSpeed );
#endif
 
 
// Handyscope HS2 Only:
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetActiveHS2) ( byte* byActiveHS2 );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetActiveHS2 GetActiveHS2;
  #endif
#else
    TP_EXTERN( word ) GetActiveHS2( byte* byActiveHS2 );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetActiveHS2) ( byte byActiveHS2 );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetActiveHS2 SetActiveHS2;
  #endif
#else
    TP_EXTERN( word ) SetActiveHS2( byte byActiveHS2 );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetDAC1451) ( double dMotor1 , double dMotor2 );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetDAC1451 SetDAC1451;
  #endif
#else
    TP_EXTERN( word ) SetDAC1451( double dMotor1 , double dMotor2 );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetMotorOn) ( byte byMotor1 , byte byMotor2 );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetMotorOn SetMotorOn;
  #endif
#else
    TP_EXTERN( word ) SetMotorOn( byte byMotor1 , byte byMotor2 );
#endif
 
 
// Obsolete routines, use is deprecated:
#ifdef TP_TYPEDEFS
    typedef dword (WINAPI * LPFN_GetMaxSampleFrequency) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetMaxSampleFrequency GetMaxSampleFrequency;
  #endif
#else
    TP_EXTERN( dword ) GetMaxSampleFrequency( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_StartMeasurement) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_StartMeasurement StartMeasurement;
  #endif
#else
    TP_EXTERN( word ) StartMeasurement( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetMeasurement) ( double* dCh1 , double* dCh2 );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetMeasurement GetMeasurement;
  #endif
#else
    TP_EXTERN( word ) GetMeasurement( double* dCh1 , double* dCh2 );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetMeasurementRaw) ( word* wCh1 , word* wCh2 );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetMeasurementRaw GetMeasurementRaw;
  #endif
#else
    TP_EXTERN( word ) GetMeasurementRaw( word* wCh1 , word* wCh2 );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetOneMeasurement) ( dword dwIndex , double* dCh1 , double* dCh2 );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetOneMeasurement GetOneMeasurement;
  #endif
#else
    TP_EXTERN( word ) GetOneMeasurement( dword dwIndex , double* dCh1 , double* dCh2 );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetOneMeasurementRaw) ( dword dwIndex , word* wCh1 , word* wCh2 );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetOneMeasurementRaw GetOneMeasurementRaw;
  #endif
#else
    TP_EXTERN( word ) GetOneMeasurementRaw( dword dwIndex , word* wCh1 , word* wCh2 );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_GetOneMeasurementRawCh) ( word wCh , dword dwIndex , word* wChData );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetOneMeasurementRawCh GetOneMeasurementRawCh;
  #endif
#else
    TP_EXTERN( word ) GetOneMeasurementRawCh( word wCh , dword dwIndex , word* wChData );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_ADC_GetData) ( word* wCh1 , word* wCh2 );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_ADC_GetData ADC_GetData;
  #endif
#else
    TP_EXTERN( word ) ADC_GetData( word* wCh1 , word* wCh2 );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_ADC_GetDataVolt) ( double* dCh1 , double* dCh2 );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_ADC_GetDataVolt ADC_GetDataVolt;
  #endif
#else
    TP_EXTERN( word ) ADC_GetDataVolt( double* dCh1 , double* dCh2 );
#endif
 
#ifdef TP_TYPEDEFS
    typedef dword (WINAPI * LPFN_GetSampleFrequency) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetSampleFrequency GetSampleFrequency;
  #endif
#else
    TP_EXTERN( dword ) GetSampleFrequency( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetSampleFrequency) ( dword* dwFreq );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetSampleFrequency SetSampleFrequency;
  #endif
#else
    TP_EXTERN( word ) SetSampleFrequency( dword* dwFreq );
#endif
 
#ifdef TP_TYPEDEFS
    typedef dword (WINAPI * LPFN_GetTriggerTimeOut) ( void );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_GetTriggerTimeOut GetTriggerTimeOut;
  #endif
#else
    TP_EXTERN( dword ) GetTriggerTimeOut( void );
#endif
 
#ifdef TP_TYPEDEFS
    typedef word (WINAPI * LPFN_SetTriggerTimeOut) ( dword dwTimeout );
  #ifndef TP_NOSTATICIMPORT
    extern LPFN_SetTriggerTimeOut SetTriggerTimeOut;
  #endif
#else
    TP_EXTERN( word ) SetTriggerTimeOut( dword dwTimeout );
#endif
 


#endif
