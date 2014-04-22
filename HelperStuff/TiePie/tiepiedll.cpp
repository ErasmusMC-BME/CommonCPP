/*
 * C file for dynamically linking a TiePie hardware dll to the application
 */


#if defined(_MSC_VER)
#include "stdafx.h"
#endif

#include "TiePieDLL.h"
#include "tiepie.h"
#include <windows.h>



//  Open / Close the instrument:
LPFN_InitInstrument InitInstrument = NULL;
LPFN_ExitInstrument ExitInstrument = NULL;
 
//  Get information about the instrument:
LPFN_GetCalibrationDate GetCalibrationDate = NULL;
LPFN_GetSerialNumber GetSerialNumber = NULL;
LPFN_GetAvailableSensitivities GetAvailableSensitivities = NULL;
LPFN_GetAvailableResolutions GetAvailableResolutions = NULL;
LPFN_GetNrChannels GetNrChannels = NULL;
LPFN_GetMaxSampleFrequencyF GetMaxSampleFrequencyF = NULL;
LPFN_GetMaxRecordLength GetMaxRecordLength = NULL;
LPFN_GetDCLevelStatus GetDCLevelStatus = NULL;
 
//  Measurement routines:
LPFN_ADC_Start ADC_Start = NULL;
LPFN_ADC_Running ADC_Running = NULL;
LPFN_ADC_Abort ADC_Abort = NULL;
LPFN_ADC_ForceTrig ADC_ForceTrig = NULL;
LPFN_ADC_Triggered ADC_Triggered = NULL;
LPFN_ADC_Ready ADC_Ready = NULL;
 
//  Retrieve the measured data:
LPFN_ADC_GetDataCh ADC_GetDataCh = NULL;
LPFN_ADC_GetDataVoltCh ADC_GetDataVoltCh = NULL;
LPFN_GetDigitalInputValues GetDigitalInputValues = NULL;
LPFN_GetOneDigitalValue GetOneDigitalValue = NULL;
 
//  Streaming measurements:
LPFN_SetDataReadyCallback SetDataReadyCallback = NULL;
LPFN_SetDataReadyEvent SetDataReadyEvent = NULL;
LPFN_SetTransferMode SetTransferMode = NULL;
LPFN_GetTransferMode GetTransferMode = NULL;
 
//  Control the input resolution in bits:
LPFN_SetResolution SetResolution = NULL;
LPFN_GetResolution GetResolution = NULL;
 
//  Control the instrument configuration:
LPFN_SetInstrumentConfig SetInstrumentConfig = NULL;
LPFN_GetInstrumentConfig GetInstrumentConfig = NULL;
 
//  Control which channels are measured:
LPFN_GetMeasureMode GetMeasureMode = NULL;
LPFN_SetMeasureMode SetMeasureMode = NULL;
 
//  Control the time base:
LPFN_GetRecordLength GetRecordLength = NULL;
LPFN_SetRecordLength SetRecordLength = NULL;
LPFN_GetPostSamples GetPostSamples = NULL;
LPFN_SetPostSamples SetPostSamples = NULL;
LPFN_GetSampleFrequencyF GetSampleFrequencyF = NULL;
LPFN_SetSampleFrequencyF SetSampleFrequencyF = NULL;
LPFN_GetExternalClock GetExternalClock = NULL;
LPFN_SetExternalClock SetExternalClock = NULL;
 
//  Control the analog input channels:
LPFN_GetSensitivity GetSensitivity = NULL;
LPFN_SetSensitivity SetSensitivity = NULL;
LPFN_GetAutoRanging GetAutoRanging = NULL;
LPFN_SetAutoRanging SetAutoRanging = NULL;
LPFN_GetCoupling GetCoupling = NULL;
LPFN_SetCoupling SetCoupling = NULL;
LPFN_GetDcLevel GetDcLevel = NULL;
LPFN_SetDcLevel SetDcLevel = NULL;
 
//  Control the trigger system:
LPFN_GetTriggerSource GetTriggerSource = NULL;
LPFN_SetTriggerSource SetTriggerSource = NULL;
LPFN_GetTriggerMode GetTriggerMode = NULL;
LPFN_SetTriggerMode SetTriggerMode = NULL;
LPFN_GetTriggerModeCh GetTriggerModeCh = NULL;
LPFN_SetTriggerModeCh SetTriggerModeCh = NULL;
LPFN_GetTriggerLevel GetTriggerLevel = NULL;
LPFN_SetTriggerLevel SetTriggerLevel = NULL;
LPFN_GetTriggerHys GetTriggerHys = NULL;
LPFN_SetTriggerHys SetTriggerHys = NULL;
LPFN_GetPXITriggerEnables GetPXITriggerEnables = NULL;
LPFN_SetPXITriggerEnables SetPXITriggerEnables = NULL;
LPFN_GetPXITriggerSlopes GetPXITriggerSlopes = NULL;
LPFN_SetPXITriggerSlopes SetPXITriggerSlopes = NULL;
 
//  Control the digital outputs:
LPFN_SetDigitalOutputs SetDigitalOutputs = NULL;
LPFN_GetDigitalOutputs GetDigitalOutputs = NULL;
 
//  Square Wave Generator:
LPFN_GetSquareWaveGenStatus GetSquareWaveGenStatus = NULL;
LPFN_GetSquareWaveGenFrequency GetSquareWaveGenFrequency = NULL;
LPFN_SetSquareWaveGenFrequency SetSquareWaveGenFrequency = NULL;
 
//  Get information about the function generator:
LPFN_GetFunctionGenStatus GetFunctionGenStatus = NULL;
LPFN_GetFuncGenMaxAmplitude GetFuncGenMaxAmplitude = NULL;
 
//  Control the Arbitrary Waveform Generator:
LPFN_GetFuncGenOutputOn GetFuncGenOutputOn = NULL;
LPFN_SetFuncGenOutputOn SetFuncGenOutputOn = NULL;
LPFN_GetFuncGenEnable GetFuncGenEnable = NULL;
LPFN_SetFuncGenEnable SetFuncGenEnable = NULL;
LPFN_GetFuncGenSignalType GetFuncGenSignalType = NULL;
LPFN_SetFuncGenSignalType SetFuncGenSignalType = NULL;
LPFN_SetFuncGenMode SetFuncGenMode = NULL;
LPFN_GetFuncGenMode GetFuncGenMode = NULL;
LPFN_GetFuncGenAmplitude GetFuncGenAmplitude = NULL;
LPFN_SetFuncGenAmplitude SetFuncGenAmplitude = NULL;
LPFN_GetFuncGenAmplitudeRange GetFuncGenAmplitudeRange = NULL;
LPFN_SetFuncGenAmplitudeRange SetFuncGenAmplitudeRange = NULL;
LPFN_GetFuncGenDCOffset GetFuncGenDCOffset = NULL;
LPFN_SetFuncGenDCOffset SetFuncGenDCOffset = NULL;
LPFN_GetFuncGenSymmetry GetFuncGenSymmetry = NULL;
LPFN_SetFuncGenSymmetry SetFuncGenSymmetry = NULL;
LPFN_GetFuncGenFrequency GetFuncGenFrequency = NULL;
LPFN_SetFuncGenFrequency SetFuncGenFrequency = NULL;
LPFN_SetFuncGenTrigSource SetFuncGenTrigSource = NULL;
LPFN_GetFuncGenTrigSource GetFuncGenTrigSource = NULL;
LPFN_FillFuncGenMemory FillFuncGenMemory = NULL;
LPFN_FuncGenBurst FuncGenBurst = NULL;
 
//  Ohm measurements ( not available by default ):
LPFN_SetupOhmMeasurements SetupOhmMeasurements = NULL;
LPFN_GetOhmValues GetOhmValues = NULL;
 
//  I2C routines:
LPFN_I2CWrite I2CWrite = NULL;
LPFN_I2CWriteNoStop I2CWriteNoStop = NULL;
LPFN_I2CRead I2CRead = NULL;
LPFN_I2CReadNoStop I2CReadNoStop = NULL;
LPFN_I2CGetSpeed I2CGetSpeed = NULL;
LPFN_I2CSetSpeed I2CSetSpeed = NULL;
 
//  Handyscope HS2 Only:
LPFN_GetActiveHS2 GetActiveHS2 = NULL;
LPFN_SetActiveHS2 SetActiveHS2 = NULL;
LPFN_SetDAC1451 SetDAC1451 = NULL;
LPFN_SetMotorOn SetMotorOn = NULL;
 
//  Obsolete routines, use is deprecated:
LPFN_GetMaxSampleFrequency GetMaxSampleFrequency = NULL;
LPFN_StartMeasurement StartMeasurement = NULL;
LPFN_GetMeasurement GetMeasurement = NULL;
LPFN_GetMeasurementRaw GetMeasurementRaw = NULL;
LPFN_GetOneMeasurement GetOneMeasurement = NULL;
LPFN_GetOneMeasurementRaw GetOneMeasurementRaw = NULL;
LPFN_GetOneMeasurementRawCh GetOneMeasurementRawCh = NULL;
LPFN_ADC_GetData ADC_GetData = NULL;
LPFN_ADC_GetDataVolt ADC_GetDataVolt = NULL;
LPFN_GetSampleFrequency GetSampleFrequency = NULL;
LPFN_SetSampleFrequency SetSampleFrequency = NULL;
LPFN_GetTriggerTimeOut GetTriggerTimeOut = NULL;
LPFN_SetTriggerTimeOut SetTriggerTimeOut = NULL;




HINSTANCE TiePieDLL;


bool OpenDLL( char *sPath ) {
  // Check is a dll was already opened. If so, close it and open the new one:
  CloseDLL();
  
  TiePieDLL = LoadLibrary( sPath );

  if (TiePieDLL != NULL) {
//  Open / Close the instrument:
    InitInstrument = (LPFN_InitInstrument) GetProcAddress( TiePieDLL , "InitInstrument" );
    ExitInstrument = (LPFN_ExitInstrument) GetProcAddress( TiePieDLL , "ExitInstrument" );
 
//  Get information about the instrument:
    GetCalibrationDate = (LPFN_GetCalibrationDate) GetProcAddress( TiePieDLL , "GetCalibrationDate" );
    GetSerialNumber = (LPFN_GetSerialNumber) GetProcAddress( TiePieDLL , "GetSerialNumber" );
    GetAvailableSensitivities = (LPFN_GetAvailableSensitivities) GetProcAddress( TiePieDLL , "GetAvailableSensitivities" );
    GetAvailableResolutions = (LPFN_GetAvailableResolutions) GetProcAddress( TiePieDLL , "GetAvailableResolutions" );
    GetNrChannels = (LPFN_GetNrChannels) GetProcAddress( TiePieDLL , "GetNrChannels" );
    GetMaxSampleFrequencyF = (LPFN_GetMaxSampleFrequencyF) GetProcAddress( TiePieDLL , "GetMaxSampleFrequencyF" );
    GetMaxRecordLength = (LPFN_GetMaxRecordLength) GetProcAddress( TiePieDLL , "GetMaxRecordLength" );
    GetDCLevelStatus = (LPFN_GetDCLevelStatus) GetProcAddress( TiePieDLL , "GetDCLevelStatus" );
 
//  Measurement routines:
    ADC_Start = (LPFN_ADC_Start) GetProcAddress( TiePieDLL , "ADC_Start" );
    ADC_Running = (LPFN_ADC_Running) GetProcAddress( TiePieDLL , "ADC_Running" );
    ADC_Abort = (LPFN_ADC_Abort) GetProcAddress( TiePieDLL , "ADC_Abort" );
    ADC_ForceTrig = (LPFN_ADC_ForceTrig) GetProcAddress( TiePieDLL , "ADC_ForceTrig" );
    ADC_Triggered = (LPFN_ADC_Triggered) GetProcAddress( TiePieDLL , "ADC_Triggered" );
    ADC_Ready = (LPFN_ADC_Ready) GetProcAddress( TiePieDLL , "ADC_Ready" );
 
//  Retrieve the measured data:
    ADC_GetDataCh = (LPFN_ADC_GetDataCh) GetProcAddress( TiePieDLL , "ADC_GetDataCh" );
    ADC_GetDataVoltCh = (LPFN_ADC_GetDataVoltCh) GetProcAddress( TiePieDLL , "ADC_GetDataVoltCh" );
    GetDigitalInputValues = (LPFN_GetDigitalInputValues) GetProcAddress( TiePieDLL , "GetDigitalInputValues" );
    GetOneDigitalValue = (LPFN_GetOneDigitalValue) GetProcAddress( TiePieDLL , "GetOneDigitalValue" );
 
//  Streaming measurements:
    SetDataReadyCallback = (LPFN_SetDataReadyCallback) GetProcAddress( TiePieDLL , "SetDataReadyCallback" );
    SetDataReadyEvent = (LPFN_SetDataReadyEvent) GetProcAddress( TiePieDLL , "SetDataReadyEvent" );
    SetTransferMode = (LPFN_SetTransferMode) GetProcAddress( TiePieDLL , "SetTransferMode" );
    GetTransferMode = (LPFN_GetTransferMode) GetProcAddress( TiePieDLL , "GetTransferMode" );
 
//  Control the input resolution in bits:
    SetResolution = (LPFN_SetResolution) GetProcAddress( TiePieDLL , "SetResolution" );
    GetResolution = (LPFN_GetResolution) GetProcAddress( TiePieDLL , "GetResolution" );
 
//  Control the instrument configuration:
    SetInstrumentConfig = (LPFN_SetInstrumentConfig) GetProcAddress( TiePieDLL , "SetInstrumentConfig" );
    GetInstrumentConfig = (LPFN_GetInstrumentConfig) GetProcAddress( TiePieDLL , "GetInstrumentConfig" );
 
//  Control which channels are measured:
    GetMeasureMode = (LPFN_GetMeasureMode) GetProcAddress( TiePieDLL , "GetMeasureMode" );
    SetMeasureMode = (LPFN_SetMeasureMode) GetProcAddress( TiePieDLL , "SetMeasureMode" );
 
//  Control the time base:
    GetRecordLength = (LPFN_GetRecordLength) GetProcAddress( TiePieDLL , "GetRecordLength" );
    SetRecordLength = (LPFN_SetRecordLength) GetProcAddress( TiePieDLL , "SetRecordLength" );
    GetPostSamples = (LPFN_GetPostSamples) GetProcAddress( TiePieDLL , "GetPostSamples" );
    SetPostSamples = (LPFN_SetPostSamples) GetProcAddress( TiePieDLL , "SetPostSamples" );
    GetSampleFrequencyF = (LPFN_GetSampleFrequencyF) GetProcAddress( TiePieDLL , "GetSampleFrequencyF" );
    SetSampleFrequencyF = (LPFN_SetSampleFrequencyF) GetProcAddress( TiePieDLL , "SetSampleFrequencyF" );
    GetExternalClock = (LPFN_GetExternalClock) GetProcAddress( TiePieDLL , "GetExternalClock" );
    SetExternalClock = (LPFN_SetExternalClock) GetProcAddress( TiePieDLL , "SetExternalClock" );
 
//  Control the analog input channels:
    GetSensitivity = (LPFN_GetSensitivity) GetProcAddress( TiePieDLL , "GetSensitivity" );
    SetSensitivity = (LPFN_SetSensitivity) GetProcAddress( TiePieDLL , "SetSensitivity" );
    GetAutoRanging = (LPFN_GetAutoRanging) GetProcAddress( TiePieDLL , "GetAutoRanging" );
    SetAutoRanging = (LPFN_SetAutoRanging) GetProcAddress( TiePieDLL , "SetAutoRanging" );
    GetCoupling = (LPFN_GetCoupling) GetProcAddress( TiePieDLL , "GetCoupling" );
    SetCoupling = (LPFN_SetCoupling) GetProcAddress( TiePieDLL , "SetCoupling" );
    GetDcLevel = (LPFN_GetDcLevel) GetProcAddress( TiePieDLL , "GetDcLevel" );
    SetDcLevel = (LPFN_SetDcLevel) GetProcAddress( TiePieDLL , "SetDcLevel" );
 
//  Control the trigger system:
    GetTriggerSource = (LPFN_GetTriggerSource) GetProcAddress( TiePieDLL , "GetTriggerSource" );
    SetTriggerSource = (LPFN_SetTriggerSource) GetProcAddress( TiePieDLL , "SetTriggerSource" );
    GetTriggerMode = (LPFN_GetTriggerMode) GetProcAddress( TiePieDLL , "GetTriggerMode" );
    SetTriggerMode = (LPFN_SetTriggerMode) GetProcAddress( TiePieDLL , "SetTriggerMode" );
    GetTriggerModeCh = (LPFN_GetTriggerModeCh) GetProcAddress( TiePieDLL , "GetTriggerModeCh" );
    SetTriggerModeCh = (LPFN_SetTriggerModeCh) GetProcAddress( TiePieDLL , "SetTriggerModeCh" );
    GetTriggerLevel = (LPFN_GetTriggerLevel) GetProcAddress( TiePieDLL , "GetTriggerLevel" );
    SetTriggerLevel = (LPFN_SetTriggerLevel) GetProcAddress( TiePieDLL , "SetTriggerLevel" );
    GetTriggerHys = (LPFN_GetTriggerHys) GetProcAddress( TiePieDLL , "GetTriggerHys" );
    SetTriggerHys = (LPFN_SetTriggerHys) GetProcAddress( TiePieDLL , "SetTriggerHys" );
    GetPXITriggerEnables = (LPFN_GetPXITriggerEnables) GetProcAddress( TiePieDLL , "GetPXITriggerEnables" );
    SetPXITriggerEnables = (LPFN_SetPXITriggerEnables) GetProcAddress( TiePieDLL , "SetPXITriggerEnables" );
    GetPXITriggerSlopes = (LPFN_GetPXITriggerSlopes) GetProcAddress( TiePieDLL , "GetPXITriggerSlopes" );
    SetPXITriggerSlopes = (LPFN_SetPXITriggerSlopes) GetProcAddress( TiePieDLL , "SetPXITriggerSlopes" );
 
//  Control the digital outputs:
    SetDigitalOutputs = (LPFN_SetDigitalOutputs) GetProcAddress( TiePieDLL , "SetDigitalOutputs" );
    GetDigitalOutputs = (LPFN_GetDigitalOutputs) GetProcAddress( TiePieDLL , "GetDigitalOutputs" );
 
//  Square Wave Generator:
    GetSquareWaveGenStatus = (LPFN_GetSquareWaveGenStatus) GetProcAddress( TiePieDLL , "GetSquareWaveGenStatus" );
    GetSquareWaveGenFrequency = (LPFN_GetSquareWaveGenFrequency) GetProcAddress( TiePieDLL , "GetSquareWaveGenFrequency" );
    SetSquareWaveGenFrequency = (LPFN_SetSquareWaveGenFrequency) GetProcAddress( TiePieDLL , "SetSquareWaveGenFrequency" );
 
//  Get information about the function generator:
    GetFunctionGenStatus = (LPFN_GetFunctionGenStatus) GetProcAddress( TiePieDLL , "GetFunctionGenStatus" );
    GetFuncGenMaxAmplitude = (LPFN_GetFuncGenMaxAmplitude) GetProcAddress( TiePieDLL , "GetFuncGenMaxAmplitude" );
 
//  Control the Arbitrary Waveform Generator:
    GetFuncGenOutputOn = (LPFN_GetFuncGenOutputOn) GetProcAddress( TiePieDLL , "GetFuncGenOutputOn" );
    SetFuncGenOutputOn = (LPFN_SetFuncGenOutputOn) GetProcAddress( TiePieDLL , "SetFuncGenOutputOn" );
    GetFuncGenEnable = (LPFN_GetFuncGenEnable) GetProcAddress( TiePieDLL , "GetFuncGenEnable" );
    SetFuncGenEnable = (LPFN_SetFuncGenEnable) GetProcAddress( TiePieDLL , "SetFuncGenEnable" );
    GetFuncGenSignalType = (LPFN_GetFuncGenSignalType) GetProcAddress( TiePieDLL , "GetFuncGenSignalType" );
    SetFuncGenSignalType = (LPFN_SetFuncGenSignalType) GetProcAddress( TiePieDLL , "SetFuncGenSignalType" );
    SetFuncGenMode = (LPFN_SetFuncGenMode) GetProcAddress( TiePieDLL , "SetFuncGenMode" );
    GetFuncGenMode = (LPFN_GetFuncGenMode) GetProcAddress( TiePieDLL , "GetFuncGenMode" );
    GetFuncGenAmplitude = (LPFN_GetFuncGenAmplitude) GetProcAddress( TiePieDLL , "GetFuncGenAmplitude" );
    SetFuncGenAmplitude = (LPFN_SetFuncGenAmplitude) GetProcAddress( TiePieDLL , "SetFuncGenAmplitude" );
    GetFuncGenAmplitudeRange = (LPFN_GetFuncGenAmplitudeRange) GetProcAddress( TiePieDLL , "GetFuncGenAmplitudeRange" );
    SetFuncGenAmplitudeRange = (LPFN_SetFuncGenAmplitudeRange) GetProcAddress( TiePieDLL , "SetFuncGenAmplitudeRange" );
    GetFuncGenDCOffset = (LPFN_GetFuncGenDCOffset) GetProcAddress( TiePieDLL , "GetFuncGenDCOffset" );
    SetFuncGenDCOffset = (LPFN_SetFuncGenDCOffset) GetProcAddress( TiePieDLL , "SetFuncGenDCOffset" );
    GetFuncGenSymmetry = (LPFN_GetFuncGenSymmetry) GetProcAddress( TiePieDLL , "GetFuncGenSymmetry" );
    SetFuncGenSymmetry = (LPFN_SetFuncGenSymmetry) GetProcAddress( TiePieDLL , "SetFuncGenSymmetry" );
    GetFuncGenFrequency = (LPFN_GetFuncGenFrequency) GetProcAddress( TiePieDLL , "GetFuncGenFrequency" );
    SetFuncGenFrequency = (LPFN_SetFuncGenFrequency) GetProcAddress( TiePieDLL , "SetFuncGenFrequency" );
    SetFuncGenTrigSource = (LPFN_SetFuncGenTrigSource) GetProcAddress( TiePieDLL , "SetFuncGenTrigSource" );
    GetFuncGenTrigSource = (LPFN_GetFuncGenTrigSource) GetProcAddress( TiePieDLL , "GetFuncGenTrigSource" );
    FillFuncGenMemory = (LPFN_FillFuncGenMemory) GetProcAddress( TiePieDLL , "FillFuncGenMemory" );
    FuncGenBurst = (LPFN_FuncGenBurst) GetProcAddress( TiePieDLL , "FuncGenBurst" );
 
//  Ohm measurements ( not available by default ):
    SetupOhmMeasurements = (LPFN_SetupOhmMeasurements) GetProcAddress( TiePieDLL , "SetupOhmMeasurements" );
    GetOhmValues = (LPFN_GetOhmValues) GetProcAddress( TiePieDLL , "GetOhmValues" );
 
//  I2C routines:
    I2CWrite = (LPFN_I2CWrite) GetProcAddress( TiePieDLL , "I2CWrite" );
    I2CWriteNoStop = (LPFN_I2CWriteNoStop) GetProcAddress( TiePieDLL , "I2CWriteNoStop" );
    I2CRead = (LPFN_I2CRead) GetProcAddress( TiePieDLL , "I2CRead" );
    I2CReadNoStop = (LPFN_I2CReadNoStop) GetProcAddress( TiePieDLL , "I2CReadNoStop" );
    I2CGetSpeed = (LPFN_I2CGetSpeed) GetProcAddress( TiePieDLL , "I2CGetSpeed" );
    I2CSetSpeed = (LPFN_I2CSetSpeed) GetProcAddress( TiePieDLL , "I2CSetSpeed" );
 
//  Handyscope HS2 Only:
    GetActiveHS2 = (LPFN_GetActiveHS2) GetProcAddress( TiePieDLL , "GetActiveHS2" );
    SetActiveHS2 = (LPFN_SetActiveHS2) GetProcAddress( TiePieDLL , "SetActiveHS2" );
    SetDAC1451 = (LPFN_SetDAC1451) GetProcAddress( TiePieDLL , "SetDAC1451" );
    SetMotorOn = (LPFN_SetMotorOn) GetProcAddress( TiePieDLL , "SetMotorOn" );
 
//  Obsolete routines, use is deprecated:
    GetMaxSampleFrequency = (LPFN_GetMaxSampleFrequency) GetProcAddress( TiePieDLL , "GetMaxSampleFrequency" );
    StartMeasurement = (LPFN_StartMeasurement) GetProcAddress( TiePieDLL , "StartMeasurement" );
    GetMeasurement = (LPFN_GetMeasurement) GetProcAddress( TiePieDLL , "GetMeasurement" );
    GetMeasurementRaw = (LPFN_GetMeasurementRaw) GetProcAddress( TiePieDLL , "GetMeasurementRaw" );
    GetOneMeasurement = (LPFN_GetOneMeasurement) GetProcAddress( TiePieDLL , "GetOneMeasurement" );
    GetOneMeasurementRaw = (LPFN_GetOneMeasurementRaw) GetProcAddress( TiePieDLL , "GetOneMeasurementRaw" );
    GetOneMeasurementRawCh = (LPFN_GetOneMeasurementRawCh) GetProcAddress( TiePieDLL , "GetOneMeasurementRawCh" );
    ADC_GetData = (LPFN_ADC_GetData) GetProcAddress( TiePieDLL , "ADC_GetData" );
    ADC_GetDataVolt = (LPFN_ADC_GetDataVolt) GetProcAddress( TiePieDLL , "ADC_GetDataVolt" );
    GetSampleFrequency = (LPFN_GetSampleFrequency) GetProcAddress( TiePieDLL , "GetSampleFrequency" );
    SetSampleFrequency = (LPFN_SetSampleFrequency) GetProcAddress( TiePieDLL , "SetSampleFrequency" );
    GetTriggerTimeOut = (LPFN_GetTriggerTimeOut) GetProcAddress( TiePieDLL , "GetTriggerTimeOut" );
    SetTriggerTimeOut = (LPFN_SetTriggerTimeOut) GetProcAddress( TiePieDLL , "SetTriggerTimeOut" );


    return( 1 );
  } else {
    return( 0 );
  } // else
}; // OpenDLL



bool CloseDLL( void ) 
{
  // Check if a dll was opened. If so, close it and open the new one:
  if (TiePieDLL != NULL) 
  {
    FreeLibrary( TiePieDLL );
    TiePieDLL = NULL;
    return( 1 );
  } // if
  else
  {
    return( 0 );
  } // else
}; // CloseDLL



