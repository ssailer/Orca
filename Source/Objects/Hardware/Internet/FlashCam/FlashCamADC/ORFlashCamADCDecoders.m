//  Orca
//  ORFlashCamADCDecoders.m
//
//  Created by Tom Caldwell on Saturday Feb 13,2021
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

#import "ORFlashCamADCDecoders.h"
#import "ORDataPacket.h"
#import "ORDataSet.h"
#import "ORDataTypeAssigner.h"
#import "ORFlashCamADCModel.h"
#import <sys/time.h>

@implementation ORFlashCamADCWaveformDecoder

- (id) init
{
    self = [super init];
    return self;
}

- (void) dealloc
{
    [fcCards release];
    [super dealloc];
}

- (void) addToObjectList:(NSMutableDictionary*)dict
{
    NSString* cname = [dict objectForKey:@"className"];
    if([cname isEqualToString:@""]) return;
    NSMutableArray* objs = [dict objectForKey:@"objects"];
    if(!objs){
        objs = [NSMutableArray array];
        [dict setObject:objs forKey:@"objects"];
    }
    [objs addObjectsFromArray:[[(ORAppDelegate*)[NSApp delegate] document]
                            collectObjectsOfClass:NSClassFromString(cname)]];
}

- (uint32_t) decodeData:(void*)someData fromDecoder:(ORDecoder*)aDecoder intoDataSet:(ORDataSet*)aDataSet
{
    /* from Model:shipEvents and Model:setWFSamples
     
     dataLengths = ((wfSamples&0xffff) << 6) | (((kFlashCamADCWFHeaderLength-3+1)&0x3f) << 22);
     dataLengths = dataLengths | ((kFlashCamADCOrcaHeaderLength&0xf) << 28);
     fprintf(stderr, "DEBUG ORFlashCamADCModel/setWFsamples: dataLength %d/0x%x\n",dataLengths, dataLengths);
     dataRecordLength = kFlashCamADCOrcaHeaderLength + (kFlashCamADCWFHeaderLength - 3 + 1) + wfSamples/2;

     
     //ship the data
     uint32_t lengths = dataLengths;
     if(!includeWF) lengths &= 0xFFFC0003F;
     fprintf(stderr, "DEBUG ORFlashCamADCModel/shipEvent: includeWF %d lengths %u/0x%x\n", includeWF, lengths, lengths);
     dataRecord[0] = dataId   | (dataRecordLength&0x3ffff);
     dataRecord[1] = lengths  | (event->type&0x3f);
     dataRecord[2] = location | ((channel&0x1f) << 9) | (index&0x1ff);
    
     */
     
    uint32* ptr = (uint32*) someData;
    uint32 length = ExtractLength(*ptr);
    
    /* from TypeAssigner.h
     
     #define kShortForm  YES
     #define kLongForm   NO

     #define kShortFormDataIdMask 0xfc000000
     #define kLongFormDataIdMask 0xfffc0000
     #define kLongFormLengthMask (~kLongFormDataIdMask)

     #define IsShortForm(x) (((x)&0x80000000) >> 31)
     #define IsLongForm(x) !IsShortForm(x)

     #define ExtractDataId(x) (IsShortForm(x) ? ((x)&kShortFormDataIdMask) : ((x)&kLongFormDataIdMask))
     #define ExtractLength(x) (IsShortForm(x) ? 1 : ((x) & ~kLongFormDataIdMask))
     */
    
    // get the header and waveform lengths, check against ExtractLength above
    uint32_t dataLengths = *(ptr+1);
    uint32_t orcaHeaderLength = (dataLengths & 0xf0000000) >> 28;
    uint32_t fcwfHeaderLength = (dataLengths & 0x0fc00000) >> 22;
    uint32_t wfSamples        = (dataLengths & 0x003fffc0) >>  6;
    // datatype is wrong here, should be int, but the info was thrown away when converting from FCIO to Orca.
    uint32_t eventType        = (dataLengths & 0x0000003f); // not used but might contain information about the type of event we have here.
    if(length != orcaHeaderLength + fcwfHeaderLength + wfSamples/2){
        NSLog(@"ORFlashCamADCWaveformDecoder: sum of orca header length %u, FCWF header length %u, and WF smaples length/2 %u != data record length %u, skipping record!\n", orcaHeaderLength, fcwfHeaderLength, wfSamples/2, length);
        return length;
    }
    
    // get the crate, card, and channel plus key strings
    uint32_t location = *(ptr+2);
    uint32_t crate   = (location & 0xf8000000) >> 27;
    uint32_t card    = (location & 0x07c00000) >> 22;
    uint32_t channel = [self getChannel:location];
    NSString* crateKey   = [self getCrateKey:crate];
    NSString* cardKey    = [self getCardKey:card];
    NSString* channelKey = [self getChannelKey:channel];
    
    // add the FPGA baseline and integrator to histograms
    ptr += orcaHeaderLength + fcwfHeaderLength - 1;
    uint16_t fpga_baseline   =  (*ptr) & 0x0000ffff;
    uint16_t fpga_integrator = ((*ptr) & 0xffff0000) >> 16;
    [aDataSet histogram:fpga_baseline numBins:0xffff sender:self
               withKeys:@"FlashCamADC", @"Baseline", crateKey, cardKey, channelKey, nil];
    [aDataSet histogram:fpga_integrator numBins:0xffff sender:self
               withKeys:@"FlashCamADC", @"Energy", crateKey, cardKey, channelKey, nil];
    
    // get the flashcam card to add to the baseline history
    NSString* key = [crateKey stringByAppendingString:cardKey];
    if(!fcCards) fcCards = [[NSMutableDictionary alloc] init];
    id obj = [fcCards objectForKey:key];
    if(!obj){
        NSMutableDictionary* dict = [NSMutableDictionary dictionaryWithObjectsAndKeys:@"ORFlashCamADCModel",
                                                                                      @"className", nil];
        [self performSelectorOnMainThread:@selector(addToObjectList:) withObject:dict waitUntilDone:YES];
        [dict setObject:@"ORFlashCamADCStdModel" forKey:@"className"];
        [self performSelectorOnMainThread:@selector(addToObjectList:) withObject:dict waitUntilDone:YES];
        NSMutableArray* listOfCards = [dict objectForKey:@"objects"];
        NSEnumerator* e = [listOfCards objectEnumerator];
        id fccard;
        while(fccard = [e nextObject]){
            if([fccard slot] == card && [fccard crateNumber] == crate){
                [fcCards setObject:fccard forKey:key];
                obj = fccard;
                break;
            }
        }
    }
    if(obj)
        if(channel>=0 && channel<[obj numberOfChannels])
            [[obj baselineHistory:channel] addDataToTimeAverage:(float)fpga_baseline];
    
    // return if the waveform is not included in this packet
    if(wfSamples == 0) return length;
    
    // only decode the waveform if it has been 100 ms since the last decoded waveform and the plotting window is open
    BOOL fullDecode = NO;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    uint64_t now = (uint64_t)(tv.tv_sec)*1000 + (uint64_t)(tv.tv_usec)/1000;
    if(!decoderOptions) decoderOptions = [[NSMutableDictionary dictionary] retain];
    NSString* lastTimeKey = [NSString stringWithFormat:@"%@,%@,%@,LastTime", crateKey, cardKey, channelKey];
    uint64_t lastTime = [[decoderOptions objectForKey:lastTimeKey] unsignedLongLongValue];
    if(now - lastTime >= 100){
        fullDecode = YES;
        [decoderOptions setObject:[NSNumber numberWithUnsignedLongLong:now] forKey:lastTimeKey];
    }
    BOOL someoneWatching = NO;
    if([aDataSet isSomeoneLooking:[NSString stringWithFormat:@"FlashCamADC,Waveforms,%d,%d,%d", crate, card, channel]]){
        someoneWatching = YES;
    }
    
    // decode the waveform if this is the first one or the above conditions are satisfied
    if(lastTime == 0 || (fullDecode && someoneWatching)){
        ptr ++;
        NSMutableData* tmpData = [NSMutableData dataWithCapacity:wfSamples*sizeof(unsigned short)];
        [tmpData setLength:wfSamples/2*sizeof(uint32_t)];
        memcpy((uint32_t*) [tmpData bytes], ptr, wfSamples*sizeof(unsigned short));
        [aDataSet loadWaveform:tmpData offset:0 unitSize:2 sender:self
                      withKeys:@"FlashCamADC", @"Waveforms", crateKey, cardKey, channelKey, nil];
    }
    
    return length;
}

- (NSString*) dataRecordDescription:(uint32_t*)dataPtr
{
    uint32_t orcaHeaderLength = (dataPtr[1] & 0xf0000000) >> 28;
    uint32_t fcwfHeaderLength = (dataPtr[1] & 0x0fc00000) >> 22;
    
    NSString* title = @"FlashCamADC Waveform Record\n\n";
    NSString* crate = [NSString stringWithFormat:@"Crate      = %U\n", (dataPtr[2] & 0xf8000000) >> 27];
    NSString* card  = [NSString stringWithFormat:@"Card       = %u\n", (dataPtr[2] & 0x07c00000) >> 22];
    NSString* chan  = [NSString stringWithFormat:@"Channel    = %u\n", [self getChannel:dataPtr[2]]];
    NSString* index = [NSString stringWithFormat:@"Ch Index   = %u\n", [self getIndex:dataPtr[2]]];
    NSString* type  = [NSString stringWithFormat:@"Event type = %u\n",  dataPtr[1] & 0x0000003f];
    uint32_t offset = orcaHeaderLength + fcwfHeaderLength - 1;
    NSString* base = [NSString stringWithFormat:@"Baseline    = %u\n",  dataPtr[offset] & 0x0000ffff];
    NSString* fint = [NSString stringWithFormat:@"Integerator = %u\n", (dataPtr[offset] & 0xffff0000) >> 16];
    offset -= fcwfHeaderLength - 1;
    NSString* header = @"Raw waveform header:\n";
    for(int i=0; i<kFlashCamADCTimeOffsetLength; i++)
        header = [header stringByAppendingFormat:@"timeoffset[%d]: %d\n", i, (int) dataPtr[offset++]];
    for(int i=0; i<kFlashCamADCDeadRegionLength; i++)
        header = [header stringByAppendingFormat:@"deadregion[%d]: %d\n", i, (int) dataPtr[offset++]];
    for(int i=0; i<kFlashCamADCTimeStampLength; i++)
        header = [header stringByAppendingFormat:@"timestamp[%d]:  %d\n", i, (int) dataPtr[offset++]];
    
    return [NSString stringWithFormat:@"%@%@%@%@%@%@%@%@%@", title, crate, card, chan, index, type, base, fint, header];
}

- (uint32_t) getChannel:(uint32_t)dataWord
{
    return (dataWord & 0x00003c00) >> 10;
}

- (uint32_t) getIndex:(uint32_t)dataWord
{
    return dataWord & 0x000003ff;
}

@end


@implementation ORFlashCamWaveformDecoder

- (id) init
{
    self = [super init];
    return self;
}

- (void) dealloc
{
    [super dealloc];
}

- (uint32_t) getChannel:(uint32_t)dataWord
{
    return (dataWord & 0x00003e00) >> 9;
}

- (uint32_t) getIndex:(uint32_t)dataWord
{
    return dataWord & 0x000001ff;
}

@end
