//  Orca
//  ORFlashCamListenerModel.h
//
//  Created by Tom Caldwell on May 1, 2020
//  Copyright (c) 2019 University of North Carolina. All rights reserved.
//-----------------------------------------------------------
//This program was prepared for the Regents of the University of
//North Carolina Department of Physics and Astrophysics
//sponsored in part by the United States
//Department of Energy (DOE) under Grant #DE-FG02-97ER41020.
//The University has certain rights in the program pursuant to
//the contract and the program should not be copied or distributed
//outside your organization.  The DOE and the University of
//North Carolina reserve all rights in the program. Neither the authors,
//University of North Carolina, or U.S. Government make any warranty,
//express or implied, or assume any liability or responsibility
//for the use of this software.
//-------------------------------------------------------------

#import "ORAuxHw.h"
#import "ORTaskSequence.h"
#import "ORDataTaker.h"
#import "ORReadOutList.h"
#import "ORTimeRate.h"
#import "ORDataFileModel.h"
#import "fcio.h"
#import "ANSIEscapeHelper.h"
#import "fsp.h"

#define kFlashCamConfigBufferLength 64
#define kFlashCamStatusBufferLength 256

@interface ORFlashCamListenerModel : ORAuxHw <ORDataTaker>
{
    @private
    NSString* interface;
    uint16_t port;
    NSString* ip;
    NSMutableArray* remoteInterfaces;
    int timeout;
    NSMutableDictionary* configParams;
    int ioBuffer;
    int stateBuffer;
    FCIOStateReader* reader;
    int readerRecordCount;
    int bufferedRecords;
    StreamProcessor* processor;
    uint32_t  configId;
    uint32_t* configBuffer;
    uint32_t  configBufferIndex;
    uint32_t  takeDataConfigIndex;
    uint32_t  bufferedConfigCount;
    uint32_t  statusId;
    uint32_t* statusBuffer;
    uint32_t  statusBufferIndex;
    uint32_t  takeDataStatusIndex;
    uint32_t  bufferedStatusCount;
    NSString* status;
    ORAlarm* runFailedAlarm;
    bool unrecognizedPacket;
    NSMutableArray* unrecognizedStates;
    NSUInteger eventCount;
    double runTime;
    double readMB;
    double rateMB;
    double rateHz;
    double timeLock;
    double deadTime;
    double totDead;
    double curDead;
    ORTimeRate* dataRateHistory;
    ORTimeRate* eventRateHistory;
    ORTimeRate* deadTimeHistory;
    ORTimeRate* swtBufferFillLevelHistory;
    ORTimeRate* swtDiscardRateHistory;
    ORTimeRate* swtOutputRateHistory;
    NSTask*     runTask;            //added. MAH 9/17/22
    NSThread* readoutThread;
    ORReadOutList* readOutList;
    NSArray* dataTakers;
    NSMutableArray* readOutArgs;
    NSMutableArray* chanMap;
    NSMutableArray* cardMap;
    bool listenerRemoteIsFile;
    int fcio_last_tag;
    bool enableStreamProcessor;
    ORDataPacket* dataPacketForThread;
    NSString* dataFileName;
    NSUInteger fclogIndex;
    NSMutableArray* fclog;
    NSMutableArray* fcrunlog;
    ORDataFileModel* dataFileObject;

    int currentStartupTime;

    NSThread* listenerThread;
    bool startFinished;
    bool setupFinished;
    bool stopRunning;
    bool takingData;
    bool runTaskCompleted;
    
    //new
    NSDateFormatter*  logDateFormatter;
    ANSIEscapeHelper* ansieHelper;

//     FSP Internal Parser variables
    int nfspHWChannels;
    int* fspHWChannelMap;
    unsigned short* fspHWPrescaleThresholds;

    int nfspPSChannels;
    int* fspPSChannelMap;
    float* fspPSChannelGains;
    float* fspPSChannelThresholds;
    int* fspPSChannelShapings;
    float* fspPSChannelLowPass;
    
    int nfspFlagChannels;
    int* fspFlagChannelMap;
    int* fspFlagChannelThresholds;
    
    int fspPulserChannel;
    int fspPulserChannelThreshold;
    int fspBaselineChannel;
    int fspBaselineChannelThreshold;
    int fspMuonChannel;
    int fspMuonChannelThreshold;
    bool debug;
}

#pragma mark •••Initialization
- (id) init;
- (id) initWithInterface:(NSString*)iface port:(uint16_t)p;
- (void) dealloc;

#pragma mark •••Accessors
- (NSString*) identifier;
- (NSString*) interface;
- (uint16_t) port;
- (NSString*) ip;
- (NSMutableArray*) remoteInterfaces;
- (NSUInteger) remoteInterfaceCount;
- (NSString*) remoteInterfaceAtIndex:(NSUInteger)index;
- (NSString*) ethType;
- (NSNumber*) configParam:(NSString*)aKey;
- (NSString*) configParamString:(NSString*)aKey;
- (uint32_t) maxADCCards;

- (NSMutableArray*) runFlags:(bool)print;
- (int) timeout;
- (int) ioBuffer;
- (int) stateBuffer;
- (FCIOStateReader*) reader;
- (int) readerRecordCount;
- (int) bufferedRecords;
- (uint32_t) configId;
- (uint32_t) statusId;
- (NSString*) status;
- (NSUInteger) eventCount;
- (double) runTime;
- (double) readMB;
- (double) rateMB;
- (double) rateHz;
- (double) timeLock;
- (double) deadTime;
- (double) totDead;
- (double) curDead;
- (double) swtRunTime;
- (int) swtEventCount;
- (double) swtAvgInputRate;
- (double) swtAvgOutputRate;
- (double) swtAvgDiscardRate;
- (int) swtFreeStates;
- (ORTimeRate*) dataRateHistory;
- (ORTimeRate*) eventRateHistory;
- (ORTimeRate*) deadTimeHistory;
- (ORTimeRate*) swtBufferFillLevelHistory;
- (ORTimeRate*) swtDiscardRateHistory;
- (ORTimeRate*) swtOutputRateHistory;
- (ORReadOutList*) readOutList;
- (NSMutableArray*) readOutArgs;
- (NSMutableArray*) children;
- (void) dataFileNameChanged:(NSNotification*) aNote;
- (NSString*) streamDescription;
- (NSUInteger) fclogLines;
- (NSString*) fclog:(NSUInteger)nprev;
- (NSUInteger) fcrunlogLines;
- (NSString*) fcrunlog:(NSUInteger)nprev;

- (void) setInterface:(NSString*)iface andPort:(uint16_t)p;
- (void) setInterface:(NSString*)iface;
- (void) updateIP;
- (void) setPort:(uint16_t)p;
- (void) setRemoteInterfaces:(NSMutableArray*)ifaces;
- (void) addRemoteInterface:(NSString*)iface;
- (void) removeRemoteInterface:(NSString*)iface;
- (void) removeRemoteInterfaceAtIndex:(NSUInteger)index;
- (void) setConfigParam:(NSString*)p withValue:(NSNumber*)v;
- (void) setConfigParam:(NSString*)p withString:(NSString*)aString;
- (void) setTimeout:(int)to;
- (void) setIObuffer:(int)io;
- (void) setStateBuffer:(int)sb;
- (void) setConfigId:(uint32_t)cId;
- (void) setStatusId:(uint32_t)sId;
- (void) setDataIds:(id)assigner;
- (void) syncDataIdsWith:(id)anotherListener;
- (void) setReadOutList:(ORReadOutList*)newList;
- (void) setReadOutArgs:(NSMutableArray*)args;
- (void) setChanMap:(NSMutableArray*)chMap;
- (void) setCardMap:(NSMutableArray*)map;
- (void) setFCLogLines:(NSUInteger)nlines;
- (void) appendToFCLog:(NSString*)line andNotify:(BOOL)notify;
- (void) clearFCLog;
- (void) appendToFCRunLog:(NSString*)line;

#pragma mark ***Formaters
- (NSDateFormatter*) logDateFormatter;
- (ANSIEscapeHelper*) ansieHelper;

#pragma mark •••Comparison methods
- (BOOL) sameInterface:(NSString*)iface andPort:(uint16_t)p;
- (BOOL) sameIP:(NSString*)address andPort:(uint16_t)p;

#pragma mark •••FCIO methods
- (bool) fcioOpen;
- (bool) fcioClose;
- (bool) fcioRead:(ORDataPacket*)aDataPacket;
- (void) runFailed;

#pragma mark •••Task methods
- (void) taskDataAvailable:(NSNotification*)note;
- (void) taskData:(NSDictionary*)taskData;
- (void) taskCompleted:(NSNotification*)note;
- (void) setupReadoutTask;
- (void) startReadoutTask;
- (void) stopReadoutTask;
- (void) readConfig:(fcio_config*)config;
- (void) readStatus:(fcio_status*)fcstatus;

#pragma mark •••Data taker methods
- (void) runTaskStarted:(ORDataPacket*)aDataPacket userInfo:(NSDictionary*)userInfo;
- (void) takeData:(ORDataPacket*)aDataPacket userInfo:(NSDictionary*)userInfo;
- (void) runIsStopping:(ORDataPacket*)aDataPacket userInfo:(NSDictionary*)userInfo;
- (void) runTaskStopped:(ORDataPacket*)aDataPacket userInfo:(NSDictionary*)userInfo;
- (void) saveReadOutList:(NSFileHandle*)aFile;
- (void) loadReadOutList:(NSFileHandle*)aFile;
- (void) reset;
- (NSDictionary*) dataRecordDescription;

#pragma mark •••Archival
- (id) initWithCoder:(NSCoder*)decoder;
- (void) encodeWithCoder:(NSCoder*)encoder;
- (NSMutableDictionary*) addParametersToDictionary:(NSMutableDictionary*)dictionary;
- (void) writeFCIOLog:(NSNotification*)note;
- (void) fileLimitExceeded:(NSNotification*)note;

- (void) runFailedMainThread;

@end

@interface ORFlashCamListenerModel (private)
- (void) setStatus:(NSString*)s;
@end

extern NSString* ORFlashCamListenerModelConfigChanged;
extern NSString* ORFlashCamListenerModelStatusChanged;
//extern NSString* ORFlashCamListenerModelConnected;
//extern NSString* ORFlashCamListenerModelDisconnected;
extern NSString* ORFlashCamListenerModelChanMapChanged;
extern NSString* ORFlashCamListenerModelCardMapChanged;
extern NSString* ORFlashCamListenerModelConfigBufferFull;
extern NSString* ORFlashCamListenerModelStatusBufferFull;
extern NSString* ORFlashCamListenerModelFCLogChanged;
extern NSString* ORFlashCamListenerModelFCRunLogChanged;
extern NSString* ORFlashCamListenerModelFCRunLogFlushed;
extern NSString* ORFlashCamListenerModelSWTConfigChanged;
extern NSString* ORFlashCamListenerModelSWTStatusChanged;
