#ifndef __APPLE__
#define _XOPEN_SOURCE 500
/* needed for random() */
#endif

#include "lpp.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <fcio.h>
#include <time_utils.h>
#include <timestamps.h>
#include <l200_channelmaps.h>


/* taken from falcon-daq*/

/*
 * Returns a pseudo-random number which is uniformly distributed in the range
 * of 0 to 1 inclusive. Uses the random(3) API, so use srandom(3) or
 * initstate(3) for seeding the PRNG.
 *
 * Note that this PRNG does not provide cryptographic-quality randomness.
 */
double random_flat(void)
{
  return random() / (double) RAND_MAX;
}


/*
 * Returns a pseudo-random number which is exponentially distributed with the
 * given mean. Uses the random(3) API, so use srandom(3) or initstate(3) for
 * seeding the PRNG.
 *
 * Note that this PRNG does not provide cryptographic-quality randomness.
 */
double random_exponential(const double mean)
{
  double u;
  while ((u = -random_flat()) == -1.0)
    ;
  return -mean * log1p(u);
}



Timestamp fcio_time_timestamps2run(int timestamps[5])
{
  assert(timestamps);

  Timestamp runtime;
  runtime.seconds = timestamps[1];
  runtime.nanoseconds = 1e9 * (double)timestamps[2] / ((double)timestamps[3] + 1.0);

  return runtime;
}

Timestamp fcio_time_ticks2run(int pps, int ticks, int maxticks)
{
  Timestamp runtime;
  runtime.seconds = pps;
  runtime.nanoseconds = 1e9 * (double)ticks / ((double)maxticks + 1.0);

  return runtime;
}

Timestamp fcio_time_run2unix(Timestamp runtime, int timeoffsets[5], int use_external_clock)
{
  assert(timeoffsets);

  if (use_external_clock) {
    runtime.seconds += timeoffsets[2];
  } else {
    runtime.seconds += timeoffsets[0];
    runtime.nanoseconds += 1000 * timeoffsets[1];
    while (runtime.nanoseconds >= 1000000000) {
      runtime.seconds++;
      runtime.nanoseconds -= 1000000000;
    }
  }
  return runtime;
}




static inline unsigned int lpp_write_decision(PostProcessor* processor, LPPState* lpp_state)
{
  if ( (lpp_state->state->last_tag != FCIOEvent) && (lpp_state->state->last_tag != FCIOSparseEvent) )
    return 1;

  if ((lpp_state->flags.event & processor->set_event_flags) || (lpp_state->flags.trigger & processor->set_trigger_flags))
    return 1;

  return 0;
}


static inline Flags lpp_evt_flags(PostProcessor* processor, Flags flags, FCIOState* state)
{

  // FCIOStateReader* reader = processor->buffer->reader;
  /*
  Determine if:
  - pulser event
  - baseline event
  - muon event
  - retrigger event
  */
  if (!state) {
    return flags;
  }

  if (state->last_tag != FCIOEvent) {
    return flags;
  }

  int nsamples = state->config->eventsamples;

  if (processor->aux.pulser_trace_index > -1) {
    if (trace_larger_than(state->event->trace[processor->aux.pulser_trace_index], 0, nsamples, nsamples, processor->aux.pulser_adc_threshold + state->event->theader[processor->aux.pulser_trace_index][0])) {
      flags.event |= EVT_AUX_PULSER;
      if (processor->loglevel >= 4) fprintf(stderr, "DEBUG lpp_evt_flags: is pulser\n");
    }
  }

  if (processor->aux.baseline_trace_index > -1) {
    if (trace_larger_than(state->event->trace[processor->aux.baseline_trace_index], 0, nsamples, nsamples, processor->aux.baseline_adc_threshold + state->event->theader[processor->aux.baseline_trace_index][0])) {
      flags.event |= EVT_AUX_BASELINE;
      if (processor->loglevel >= 4) fprintf(stderr, "DEBUG lpp_evt_flags: is baseline\n");
    }
  }

  if (processor->aux.muon_trace_index > -1) {
    if (trace_larger_than(state->event->trace[processor->aux.muon_trace_index], 0, nsamples, nsamples, processor->aux.muon_adc_threshold + state->event->theader[processor->aux.muon_trace_index][0])) {
      flags.event |= EVT_AUX_MUON;
      if (processor->loglevel >= 4) fprintf(stderr, "DEBUG lpp_evt_flags: is muon\n");
    }
  }

  LPPState* previous_lpp_state = NULL;
  int previous_counter = 0;
  while ( (previous_lpp_state = LPPBufferGetState(processor->buffer,previous_counter--)) ) {
    if (!previous_lpp_state || !previous_lpp_state->in_buffer)
      break;

    FCIOState* previous_state = previous_lpp_state->state;
    /* TODO Check if Sparseevents retrigger, and how to handle this case*/
    if (previous_state->last_tag != FCIOEvent)
      continue;

    Timestamp now_ts = fcio_time_timestamps2run(state->event->timestamp);
    Timestamp previous_ts = fcio_time_timestamps2run(previous_state->event->timestamp);

    Timestamp delta_ts = { .seconds = 0, .nanoseconds = 0 };
    if (state->config->adcbits == 16) {
      delta_ts.seconds = 0;
      delta_ts.nanoseconds = state->config->eventsamples * 16;
    } else if (state->config->adcbits == 12){
      delta_ts.seconds = 0;
      delta_ts.nanoseconds = state->config->eventsamples * 4;
    } else {
      fprintf(stderr, "CRITICAL lpp_evt_flags: Only support 12- or 16-bit ADC data. Got %d bit precision.\n", state->config->adcbits);
    }
    Timestamp event_delta = timestamp_sub(now_ts, previous_ts);
    if (timestamp_leq(event_delta, delta_ts)) {
      flags.event |= EVT_RETRIGGER;
      /* Marking the previous event as EVT_EXTENDED happens in UpdateLPPState. */
      if (processor->loglevel >= 4) fprintf(stderr, "DEBUG lpp_evt_flags: is retrigger with now=%ld.%09ld previous=%ld.%09ld delta=%ld.%09ld\n",
        now_ts.seconds, now_ts.nanoseconds,
        previous_ts.seconds, previous_ts.nanoseconds,
        delta_ts.seconds, delta_ts.nanoseconds
      );
    }

  }

  return flags;
}


static inline Flags lpp_st_majority(PostProcessor* processor, Flags flags, FCIOState* state)
{
  fcio_config* config = state->config;
  fcio_event* event = state->event;

  tale_dsp_fpga_energy_majority(processor->fpga_majority_cfg, config->adcs, event->theader);

  /* don't set to higher than 1 if you are sane */
  if (processor->fpga_majority_cfg->majority >= processor->majority_threshold) {
    flags.event |= EVT_FPGA_MULTIPLICITY;

    /* if majority is >= 1, then the following check is safe, otherwise think about what happens when majority == 0
       if there is any channel with a majority above the threshold, it's a force trigger, if not, it should be
       prescaled and not affect the rest of the datastream.
    */
    if (processor->fpga_majority_cfg->n_fpga_energy_below == processor->fpga_majority_cfg->majority)
      flags.event |= EVT_FPGA_MULTIPLICITY_ENERGY_BELOW;
    else
      flags.trigger |= ST_TRIGGER_FORCE;
  }

  return flags;
}

static inline Flags lpp_st_analogue_sum(PostProcessor* processor, Flags flags, FCIOState* state)
{
  fcio_config* config = state->config;
  fcio_event* event = state->event;

  tale_dsp_windowed_analogue_sum(processor->analogue_sum_cfg, config->eventsamples, config->adcs, event->trace);

  if (processor->analogue_sum_cfg->largest_sum_pe >= processor->sum_threshold_pe) {
    flags.trigger |= ST_TRIGGER_SIPM_NPE;
  }

  if (processor->analogue_sum_cfg->largest_sum_pe >= processor->windowed_sum_threshold_pe) {
    flags.event |= EVT_ASUM_MIN_NPE;
  }

  if (processor->muon_coincidence && flags.event & EVT_AUX_MUON) {
    flags.trigger |= ST_TRIGGER_FORCE;
  }

  return flags;
}

static inline Flags lpp_st_prescaling(PostProcessor* processor, Flags flags, Timestamp event_unix_timestamp)
{
  /*
    Exit early if:
    - prescaling has not been turned on
    - the event categority contains non-phy triggers / occasions
    - if the software trigger found a force trigger

    -> should result in pure veto triggers being prescaled.
  */

  /* ge prescaling */
  if ( processor->ge_prescaling_rate > 0.0
       && (flags.event & EVT_FPGA_MULTIPLICITY_ENERGY_BELOW)
       && ((flags.trigger & ST_TRIGGER_FORCE) == 0)
    ) {
    if (processor->ge_prescaling_timestamp.seconds == 0) {
      /* initialize with the first event in the stream, which will have the prescaled flag.*/
      processor->ge_prescaling_timestamp = event_unix_timestamp;
    }
    if (timestamp_geq(event_unix_timestamp, processor->ge_prescaling_timestamp) ) {
      flags.trigger |= ST_TRIGGER_GE_PRESCALED;
      double shift = random_exponential(1.0 / processor->ge_prescaling_rate);
      double integral;
      double fractional = modf(shift, &integral);
      processor->ge_prescaling_timestamp.seconds += (long)integral;
      processor->ge_prescaling_timestamp.nanoseconds += (long)(fractional*1.0e9);
    }
  }

  /* sipm prescaling
    Only prescales events which have not otherwise been triggered.
    Ge prescaling takes precedence.
    TODO: Should this be a combined rate, i.e. the sipm_prescaling_timestamp be shifted even if it was 
          a ge_prescaled event?
  */
  if ( !processor->sipm_prescaling || flags.trigger )
    return flags;

  switch (processor->sipm_prescaling[0]) {
    case 'A': {
      if (processor->sipm_prescaling_timestamp.seconds == 0) {
        /* initialize with the first event in the stream, which will have the prescaled flag.*/
        processor->sipm_prescaling_timestamp = event_unix_timestamp;
      }

      if (timestamp_geq(event_unix_timestamp, processor->sipm_prescaling_timestamp) ) {
        flags.trigger |= ST_TRIGGER_SIPM_PRESCALED;
        double shift = random_exponential(1.0 / processor->sipm_prescaling_rate);
        double integral;
        double fractional = modf(shift, &integral);
        processor->sipm_prescaling_timestamp.seconds += (long)integral;
        processor->sipm_prescaling_timestamp.nanoseconds += (long)(fractional*1.0e9);
      }

      break;
    }
    case 'F': {
      if (processor->sipm_prescaling_counter == processor->sipm_prescaling_offset) {
        flags.trigger |= ST_TRIGGER_SIPM_PRESCALED;
        processor->sipm_prescaling_counter = 0;

      } else {
        processor->sipm_prescaling_counter++;

      }

      break;
    }
  }
  return flags;
}



int lpp_process_fcio_state(PostProcessor* processor, LPPState* lpp_state, FCIOState* state)
{

  AnalogueSumCfg* analogue_sum_cfg = processor->analogue_sum_cfg;
  FPGAMajorityCfg* fpga_majority_cfg = processor->fpga_majority_cfg;

  lpp_state->stream_tag = state->last_tag;

  switch (state->last_tag)
  {
    case FCIOEvent:
    case FCIOSparseEvent: {

      lpp_state->timestamp = fcio_time_timestamps2run(state->event->timestamp);
      lpp_state->unixstamp = fcio_time_run2unix(lpp_state->timestamp, state->event->timeoffset, state->config->gps);

      Flags flags = { .event = EVT_NULL, .trigger = ST_NULL};

      if (processor->checks) {

        const int max_ticks = 249999999;
        if (state->event->timestamp[2] > max_ticks) {
          if (processor->loglevel >= 2) fprintf(stderr, "WARNING timestamp of event %i is out of bounds (ticks=%u)\n", state->event->timestamp[0], state->event->timestamp[2]);
        }

        if (state->event->timestamp[3] != max_ticks) {
          if (processor->loglevel >= 2) fprintf(stderr, "WARNING lost time synchronisation in previous second (max_ticks=%d)\n", state->event->timestamp[3]);
        }
      }
      flags = lpp_evt_flags(processor, flags, state);

      if (fpga_majority_cfg) {
        flags = lpp_st_majority(processor, flags, state);
        lpp_state->majority = fpga_majority_cfg->majority;
        lpp_state->ge_max_fpga_energy = fpga_majority_cfg->max_fpga_energy;
        lpp_state->ge_min_fpga_energy = fpga_majority_cfg->min_fpga_energy;
      }

      if (analogue_sum_cfg) {
        flags = lpp_st_analogue_sum(processor, flags, state);
        lpp_state->largest_sum_pe = analogue_sum_cfg->largest_sum_pe;
        lpp_state->largest_sum_offset = analogue_sum_cfg->largest_sum_offset;
        lpp_state->largest_pe = analogue_sum_cfg->largest_pe;
        lpp_state->channel_multiplicity = analogue_sum_cfg->multiplicity;
      }

      flags = lpp_st_prescaling(processor, flags, lpp_state->unixstamp);

      lpp_state->flags = flags;

      lpp_state->contains_timestamp = 1;
      break;
    }

    case FCIORecEvent: {
      lpp_state->timestamp = fcio_time_timestamps2run(state->recevent->timestamp);
      lpp_state->unixstamp = fcio_time_run2unix(lpp_state->timestamp, state->event->timeoffset, state->config->gps);

      lpp_state->contains_timestamp = 1;

      break;
    }



    case FCIOConfig: {
      // format == 0 is already converted to trace indices used by fcio
      if (processor->aux.tracemap_format) {
        convert_trace_idx(1, &processor->aux.pulser_trace_index, processor->aux.tracemap_format, state->config->tracemap);
        convert_trace_idx(1, &processor->aux.baseline_trace_index, processor->aux.tracemap_format, state->config->tracemap);
        convert_trace_idx(1, &processor->aux.muon_trace_index, processor->aux.tracemap_format, state->config->tracemap);
      }

      if (processor->analogue_sum_cfg) {
        if (processor->analogue_sum_cfg->tracemap_format) {
          convert_trace_idx(processor->analogue_sum_cfg->ntraces, processor->analogue_sum_cfg->tracemap, processor->analogue_sum_cfg->tracemap_format, state->config->tracemap);
        }

        if (analogue_sum_cfg->sum_window_stop_sample < 0)
          analogue_sum_cfg->sum_window_stop_sample = state->config->eventsamples;

        if (analogue_sum_cfg->sum_window_stop_sample + analogue_sum_cfg->dsp_post_max_samples > state->config->eventsamples) {
          analogue_sum_cfg->sum_window_stop_sample = state->config->eventsamples - analogue_sum_cfg->dsp_post_max_samples;
        }
        if (analogue_sum_cfg->sum_window_start_sample - analogue_sum_cfg->dsp_pre_max_samples < 0) {
          analogue_sum_cfg->sum_window_start_sample = analogue_sum_cfg->dsp_pre_max_samples;
        }
        int valid_window = analogue_sum_cfg->sum_window_stop_sample - analogue_sum_cfg->sum_window_start_sample;
        if (valid_window <= 0) {
          fprintf(stderr, "CRITICAL lpp_process_fcio_state: sum_window_start_sample %d and sum_window_stop_sample %d overlap! No samples will be checked.\n", analogue_sum_cfg->sum_window_start_sample, analogue_sum_cfg->sum_window_stop_sample);
        }
        if (valid_window < analogue_sum_cfg->coincidence_window) {
          if (processor->loglevel)
            fprintf(stderr, "ERROR lpp_process_fcio_state: not enough samples for these dsp settings to allow a analog sum window of size %d samples, reduced to %d.\n", analogue_sum_cfg->coincidence_window, valid_window);
          analogue_sum_cfg->coincidence_window = valid_window;
        }

        for (int i = 0; i < analogue_sum_cfg->ntraces; i++) {
          analogue_sum_cfg->dsp_start_sample[i] = analogue_sum_cfg->sum_window_start_sample - analogue_sum_cfg->dsp_pre_samples[i];
          analogue_sum_cfg->dsp_stop_sample[i] = analogue_sum_cfg->sum_window_stop_sample + analogue_sum_cfg->dsp_post_samples[i];
          if (processor->loglevel >= 4) {
            fprintf(stderr, "DEBUG lpp_process_fcio_state: adjusting windows: channel %d sum start %d sum stop %d dsp start %d dsp stop %d\n", i,
              analogue_sum_cfg->sum_window_start_sample, analogue_sum_cfg->sum_window_stop_sample,
              analogue_sum_cfg->dsp_start_sample[i], analogue_sum_cfg->dsp_stop_sample[i]
            );
          }
        }
      }

      if (processor->fpga_majority_cfg) {
        if (processor->fpga_majority_cfg->tracemap_format) {
          convert_trace_idx(processor->fpga_majority_cfg->ntraces, processor->fpga_majority_cfg->tracemap, processor->fpga_majority_cfg->tracemap_format, state->config->tracemap);
        }

        if (fpga_majority_cfg->ntraces <= 0) {
          fpga_majority_cfg->ntraces = state->config->adcs;
          processor->majority_threshold = 0; // we pass the event even if there is no fpga_energy; it must be in the stream for a reason.
        }

      }
      lpp_state->contains_timestamp = 0;

      break;
    }
    case FCIOStatus:
      // TODO check if the the data[0].pps/ticks/maxticks is a valid timestamp for the statuspaket
      // the status.statustime throws away the nanoseconds. The only advantage would be, that we might
      // send some soft triggered events a bit earlier, as the pre-trigger timestamp shift, if we use the
      // status timestamp here.
      // lpp_state->timestamp = fcio_time_ticks2run(lpp_state->state->status->statustime); break;
    default: {
      lpp_state->contains_timestamp = 0;
    break;
    }
  }

  return lpp_state->contains_timestamp;
}

void lpp_process_timings(PostProcessor* processor, LPPState* lpp_state)
{

  if (timestamp_geq(processor->post_trigger_timestamp,lpp_state->timestamp)) {

    lpp_state->flags.event |= EVT_FORCE_POST_WINDOW;

    if (lpp_state->flags.event & EVT_ASUM_MIN_NPE) {
      lpp_state->flags.trigger |= ST_TRIGGER_SIPM_NPE_IN_WINDOW;
    }
  }

  if (lpp_state->flags.trigger & ST_TRIGGER_FORCE) {
    /* current state is germanium trigger, keep it and start checking all previous and future states against
      the trigger windows.
    */
    processor->force_trigger_timestamp = lpp_state->timestamp;
    processor->post_trigger_timestamp = timestamp_add(processor->force_trigger_timestamp, processor->post_trigger_window);
    processor->pre_trigger_timestamp = timestamp_sub(processor->force_trigger_timestamp, processor->pre_trigger_window);

    LPPState* update_lpp_state = NULL;
    int previous_counter = 0; // current lpp_state is a peeked state, so GetState(0) is the "previous one"
    while ( (update_lpp_state = LPPBufferGetState(processor->buffer, previous_counter--))
            && timestamp_geq(update_lpp_state->timestamp, processor->pre_trigger_timestamp)
          ) {
            update_lpp_state->flags.event |= EVT_FORCE_PRE_WINDOW;
            if (update_lpp_state->flags.event & EVT_ASUM_MIN_NPE) {
              update_lpp_state->flags.trigger |= ST_TRIGGER_SIPM_NPE_IN_WINDOW;
            }
    }
  }


  if (lpp_state->flags.event & EVT_RETRIGGER) {
    LPPState* previous_lpp_state = LPPBufferGetState(processor->buffer, 0);
    // fprintf(stderr, "DEBUDEBUG %ld %d %d\n", lpp_state - previous_lpp_state, lpp_state->state->event->timestamp[0], previous_lpp_state->state->event->timestamp[0]);
    if (previous_lpp_state && ((previous_lpp_state->flags.event & EVT_RETRIGGER) == 0)) {
      // fprintf(stderr, "setting %d to E\n", previous_lpp_state->state->event->timestamp[0]);
      previous_lpp_state->flags.event |= EVT_EXTENDED;
    }
  }

}


void lpp_process(PostProcessor* processor, LPPState* lpp_state, FCIOState* state)
{
  lpp_state->state = state;

  if (lpp_process_fcio_state(processor, lpp_state, state)) {
    lpp_process_timings(processor, lpp_state);

  }

}

int LPPInput(PostProcessor* processor, FCIOState* state)
{
  if (!processor || !state)
    return 0;

  LPPState* lpp_state = LPPBufferPeek(processor->buffer);

  if (!lpp_state) {
    fprintf(stderr, "CRITICAL LPPInput: Buffer full, refuse to overwrite.\n");
    return 0;
  }
  if (state->last_tag == FCIOEvent || state->last_tag == FCIOSparseEvent ) {
    processor->nevents_read++;
  }
  processor->nrecords_read++;

  lpp_process(processor, lpp_state, state);

  LPPBufferCommit(processor->buffer);

  return processor->buffer->max_states - processor->buffer->fill_level;
  /* could return possible number of state to write here, allowing know how often to call the LPPGetNextState, would be also the used/needed buffer depth*/
}

LPPState* LPPOutput(PostProcessor* processor)
{
  if (!processor)
    return NULL;

  LPPState* lpp_state = LPPBufferFetch(processor->buffer);

  if (!lpp_state) {
    return NULL;
  }

  lpp_state->write = lpp_write_decision(processor, lpp_state);

  if (lpp_state->write) {
    processor->nrecords_written++;
    if ((lpp_state->state->last_tag == FCIOSparseEvent) || (lpp_state->state->last_tag == FCIOEvent))
      processor->nevents_written++;

  } else {
    processor->nrecords_discarded++;
    if ((lpp_state->state->last_tag == FCIOSparseEvent) || (lpp_state->state->last_tag == FCIOEvent))
      processor->nevents_discarded++;

  }



  return lpp_state;

}

void LPPEnableTriggerFlags(PostProcessor* processor, unsigned int flags)
{
  processor->set_trigger_flags = flags;
  if (processor->loglevel >= 4)
    fprintf(stderr, "DEBUG LPPEnableTriggerFlags: %u\n", flags);
}

void LPPEnableEventFlags(PostProcessor* processor, unsigned int flags)
{
  processor->set_event_flags = flags;
  if (processor->loglevel >= 4)
    fprintf(stderr, "DEBUG LPPEnableEventFlags: %u\n", flags);
}

PostProcessor* LPPCreate()
{
  PostProcessor* processor = calloc(1, sizeof(PostProcessor));

  processor->minimum_buffer_window.seconds = 0;
  processor->minimum_buffer_window.nanoseconds = (FCIOMaxSamples+1) * 16; // this is required to check for retrigger events
  processor->minimum_buffer_depth = 2; // we accept one FCIOStatus packet interleaved between two events
  processor->stats.start_time = 0.0; // reset, actual start time happens with the first record insertion.

  /* hardcoded defaults which should make sense. Used SetFunctions outside to overwrite */
  LPPEnableEventFlags(processor, EVT_AUX_PULSER | EVT_AUX_BASELINE | EVT_EXTENDED | EVT_RETRIGGER );
  LPPEnableTriggerFlags(processor, ST_TRIGGER_FORCE | ST_TRIGGER_SIPM_NPE_IN_WINDOW | ST_TRIGGER_SIPM_NPE | ST_TRIGGER_SIPM_PRESCALED | ST_TRIGGER_GE_PRESCALED );

  return processor;
}

void LPPSetLogLevel(PostProcessor* processor, int loglevel)
{
  processor->loglevel = loglevel;
}

void LPPSetLogTime(PostProcessor* processor, double log_time)
{
  processor->stats.log_time = log_time;
}


int LPPSetBufferSize(PostProcessor* processor, int buffer_depth)
{
  if (processor->buffer) {
    LPPBufferDestroy(processor->buffer);
  }
  if (buffer_depth < processor->minimum_buffer_depth)
    buffer_depth = processor->minimum_buffer_depth;


  Timestamp buffer_window = timestamp_greater(processor->minimum_buffer_window, processor->pre_trigger_window) ? processor->minimum_buffer_window : processor->pre_trigger_window;
  processor->buffer = LPPBufferCreate(buffer_depth, buffer_window);
  if (!processor->buffer) {
    if (processor->loglevel)
      fprintf(stderr, "ERROR LPPSetBufferSize: Couldn't allocate LPPBuffer.\n");

    return 0;
  }
  return 1;
}

void LPPDestroy(PostProcessor* processor)
{
  LPPBufferDestroy(processor->buffer);
  free(processor->fpga_majority_cfg);
  free(processor->analogue_sum_cfg);
  free(processor);
}


int LPPFlush(PostProcessor* processor)
{
  if (!processor)
    return 0;

  return LPPBufferFlush(processor->buffer);
}

int LPPFreeStates(PostProcessor* processor)
{
  if (!processor)
    return 0;

  return processor->buffer->max_states - processor->buffer->fill_level;
}

static inline void lpp_init_stats(PostProcessor* processor) {
  LPPStats* stats = &processor->stats;
  if (stats->start_time == 0.0) {
    stats->dt_logtime = stats->start_time = elapsed_time(0.0);
  }
}


int LPPStatsUpdate(PostProcessor* processor, int force)
{
  LPPStats* stats = &processor->stats;

  if ( elapsed_time(stats->dt_logtime) > stats->log_time || force) {
    stats->runtime = elapsed_time(stats->start_time);
    stats->dt_current = elapsed_time(stats->dt_logtime);
    stats->dt_logtime = elapsed_time(0.0);

    stats->current_nread = processor->nevents_read - stats->n_read_events;
    stats->n_read_events = processor->nevents_read;

    stats->current_nwritten = processor->nevents_written - stats->n_written_events;
    stats->n_written_events = processor->nevents_written;

    stats->n_discarded_events = stats->n_read_events - stats->n_written_events;

    stats->current_read_rate = stats->current_nread / stats->dt_current;
    stats->current_write_rate = stats->current_nwritten / stats->dt_current;
    stats->current_discard_rate = stats->current_read_rate - stats->current_write_rate;

    stats->avg_read_rate = processor->nevents_read / stats->runtime;
    stats->avg_write_rate = processor->nevents_written / stats->runtime;
    stats->avg_discard_rate = stats->avg_read_rate - stats->avg_write_rate;

    return 1;
  }
  return 0;
}

int LPPStatsInfluxString(PostProcessor* processor, char* logstring, size_t logstring_size)
{

//  if (LPPStatsUpdate(processor, !lpp_state)) {

    LPPStats* stats = &processor->stats;

    int ret = snprintf(logstring, logstring_size,"run_time=%.03f,cur_read=%.03f,cur_write=%.03f,cur_discard=%.03f,avg_read=%.03f,avg_write=%.03f,avg_discard=%.03f,cur_nread=%d,cur_nsent=%d,total_nread=%d,total_nsent=%d,total_ndiscarded=%d",
      stats->runtime,
      stats->current_read_rate, stats->current_write_rate, stats->current_discard_rate,
      stats->avg_read_rate, stats->avg_write_rate, stats->avg_discard_rate,
      stats->current_nread, stats->current_nwritten,
      stats->n_read_events,
      stats->n_written_events,
      stats->n_discarded_events
    );

  if (ret >= 0 && ret < (int)logstring_size)
    return 1;
//  }
  return 0;
}

LPPState* LPPGetNextState(PostProcessor* processor, FCIOStateReader* reader, int* timedout)
{
  /*
  - check if buffered -> write
  - no buffered, getnextstate until buffered reader -> write
  - no buffered and getnextstate null -> exit
  */
  if (timedout)
    *timedout = 0;

  if (!processor || !reader)
    return NULL;

  if (!processor->buffer && LPPSetBufferSize(processor, reader->max_states-1))
    return NULL;

  lpp_init_stats(processor);


  LPPState* lpp_state;
  FCIOState* state;
  int nfree = LPPFreeStates(processor);

  while ( !(lpp_state = LPPOutput(processor)) ) {
    /* this should handle the timeout, so we don't have to do it in the postprocessor.
     */
    if (!nfree) {
      if (timedout)
        *timedout = 3;
      return NULL;
    }
    state = FCIOGetNextState(reader, timedout);

    if (!state ) {
      /* End-of-Stream or timeout reached, we assume finish and flush
        timedout == 0 here means stream is closed, no timeout was reached
        timedout == 1 no event from buffer for specified timeout interval.
        timedout == 2 stream is still alive, but only unselected tags arrived
      */
      if (LPPFlush(processor)) {
        continue; // still something in the buffer, try to read the rest after unlocking the buffer window.
      } else {
        return NULL;
      }
    } else {
      // int n_free_buffer_states = LPPInput(processor, state); // got a valid state, process and hope that we get a new lpp_state
      nfree = LPPInput(processor, state); // got a valid state, process and hope that we get a new lpp_state on the next loop
    }
  }

  return lpp_state;
}





int LPPSetAuxParameters(PostProcessor* processor, const char* channelmap_format,
                                              int digital_pulser_channel, int pulser_level_adc,
                                              int digital_baseline_channel, int baseline_level_adc,
                                              int digital_muon_channel, int muon_level_adc
                                              )
{
  if ((processor->aux.tracemap_format = get_channelmap_format(channelmap_format)) < 0) {
    if (processor->loglevel)
      fprintf(stderr, "ERROR LPPSetAuxParameters: channel map type %s is not supported. Valid inputs are \"fcio-trace-index\", \"fcio-tracemap\" or \"rawid\".\n", channelmap_format);
    return 0;
  }

  processor->aux.pulser_trace_index = digital_pulser_channel;
  processor->aux.pulser_adc_threshold = pulser_level_adc;
  processor->aux.baseline_trace_index = digital_baseline_channel;
  processor->aux.baseline_adc_threshold = baseline_level_adc;
  processor->aux.muon_trace_index = digital_muon_channel;
  processor->aux.muon_adc_threshold = muon_level_adc;


  if (processor->loglevel >= 4) {
    fprintf(stderr, "DEBUG LPPSetAuxParameters\n");
    fprintf(stderr, "DEBUG digital_pulser_channel   %d level_adc %d\n", processor->aux.pulser_trace_index, processor->aux.pulser_adc_threshold);
    fprintf(stderr, "DEBUG digital_baseline_channel %d level_adc %d\n", processor->aux.baseline_trace_index, processor->aux.baseline_adc_threshold);
    fprintf(stderr, "DEBUG digital_muon_channel     %d level_adc %d\n", processor->aux.muon_trace_index, processor->aux.muon_adc_threshold);

  }
  return 1;
}

int LPPSetGeParameters(PostProcessor* processor, int nchannels, int* channelmap, const char* channelmap_format, int majority_threshold, int skip_full_counting, unsigned short* ge_prescaling_threshold_adc, float ge_average_prescaling_rate_hz)
{
  processor->fpga_majority_cfg = calloc(1, sizeof(FPGAMajorityCfg));
  FPGAMajorityCfg* fmc = processor->fpga_majority_cfg;

  if ((fmc->tracemap_format = get_channelmap_format(channelmap_format)) < 0) {
    if (processor->loglevel)
      fprintf(stderr, "ERROR LPPSetGeParameters: channel map type %s is not supported. Valid inputs are \"fcio-trace-index\", \"fcio-tracemap\" or \"rawid\".\n", channelmap_format);
    free(fmc);
    return 0;
  }
  fmc->ntraces = nchannels;

  for (int i = 0; i < nchannels && i < FCIOMaxChannels; i++) {
    fmc->tracemap[i] = channelmap[i];
    fmc->fpga_energy_threshold_adc[i] = ge_prescaling_threshold_adc[i];
  }
  fmc->fast = skip_full_counting;
  processor->majority_threshold = majority_threshold;
  processor->ge_prescaling_rate = ge_average_prescaling_rate_hz;

  if (processor->loglevel >= 4) {
    fprintf(stderr, "DEBUG LPPSetGeParameters\n");
    fprintf(stderr, "DEBUG majority_threshold %d\n", processor->majority_threshold);
    fprintf(stderr, "DEBUG skip_full_counting %d\n", fmc->fast);
    for (int i = 0; i < fmc->ntraces; i++)
      fprintf(stderr, "DEBUG channel %d\n", fmc->tracemap[i]);
  }
  return 1;
}



int LPPSetSiPMParameters(PostProcessor* processor, int nchannels, int* channelmap, const char* channelmap_format,
                                                float* calibration_pe_adc, float* channel_thresholds_pe,
                                                int* shaping_width_samples, float* lowpass_factors,
                                                int coincidence_pre_window_ns, int coincidence_post_window_ns,
                                                int coincidence_window_samples,
                                                int sum_window_start_sample, int sum_window_stop_sample,
                                                float sum_threshold_pe, float coincidence_sum_threshold_pe,
                                                float average_prescaling_rate_hz, int enable_muon_coincidence
                                                )
{
  processor->analogue_sum_cfg = calloc(1, sizeof(AnalogueSumCfg));
  AnalogueSumCfg* asc = processor->analogue_sum_cfg;

  if ((asc->tracemap_format = get_channelmap_format(channelmap_format)) < 0) {
    if (processor->loglevel)
      fprintf(stderr, "ERROR LPPSetSiPMParameters: channel map type %s is not supported. Valid inputs are \"fcio-trace-index\", \"fcio-tracemap\" or \"rawid\".\n", channelmap_format);
    free(asc);
    return 0;
  }

  processor->windowed_sum_threshold_pe = coincidence_sum_threshold_pe;
  processor->sum_threshold_pe = sum_threshold_pe;
  processor->pre_trigger_window.seconds = 0;
  processor->pre_trigger_window.nanoseconds = coincidence_pre_window_ns;
  processor->post_trigger_window.seconds = 0;
  processor->post_trigger_window.nanoseconds = coincidence_post_window_ns;
  processor->sipm_prescaling_rate = average_prescaling_rate_hz;
  if (processor->sipm_prescaling_rate == 0.0)
    processor->sipm_prescaling = "AverageRate"; // could be "FixedOffset" when selecting ->sipm_prescaling_offset, but is disabled here.
  else
    processor->sipm_prescaling = NULL;
  processor->muon_coincidence = enable_muon_coincidence;

  // asc->repetition = processor->fast?4:sma_repetition;
  asc->coincidence_window = coincidence_window_samples;
  asc->sum_window_start_sample = sum_window_start_sample;
  asc->sum_window_stop_sample = sum_window_stop_sample;

  /* TODO CHECK THIS*/
  asc->apply_gain_scaling = 1;


  asc->ntraces = nchannels;

  asc->dsp_pre_max_samples = 0;
  asc->dsp_post_max_samples = 0;
  for (int i = 0; i < nchannels && i < FCIOMaxChannels; i++) {
    asc->tracemap[i] = channelmap[i];
    asc->gains[i] = calibration_pe_adc[i];
    asc->thresholds[i] = channel_thresholds_pe[i];
    asc->shaping_widths[i] = shaping_width_samples[i];
    asc->lowpass[i] = lowpass_factors[i];
    asc->dsp_pre_samples[i] = tale_dsp_diff_and_smooth_pre_samples(shaping_width_samples[i], asc->lowpass[i]);
    if (asc->dsp_pre_samples[i] > asc->dsp_pre_max_samples)
      asc->dsp_pre_max_samples = asc->dsp_pre_samples[i];
    asc->dsp_post_samples[i] = tale_dsp_diff_and_smooth_post_samples(shaping_width_samples[i], asc->lowpass[i]);
    if (asc->dsp_post_samples[i] > asc->dsp_post_max_samples)
      asc->dsp_post_max_samples = asc->dsp_post_samples[i];
  }

  if (processor->loglevel >= 4) {
    /* DEBUGGING enabled, print all inputs */
    fprintf(stderr, "DEBUG LPPSetSiPMParameters:\n");
    fprintf(stderr, "DEBUG average_prescaling_rate_hz   %f\n", processor->sipm_prescaling_rate);
    fprintf(stderr, "DEBUG sum_window_start_sample      %d\n", asc->sum_window_start_sample);
    fprintf(stderr, "DEBUG sum_window_stop_sample       %d\n", asc->sum_window_stop_sample);
    fprintf(stderr, "DEBUG dsp_pre_max_samples          %d\n", asc->dsp_pre_max_samples);
    fprintf(stderr, "DEBUG dsp_post_max_samples         %d\n", asc->dsp_post_max_samples);
    fprintf(stderr, "DEBUG coincidence_pre_window_ns    %ld\n", processor->pre_trigger_window.nanoseconds);
    fprintf(stderr, "DEBUG coincidence_post_window_ns   %ld\n", processor->post_trigger_window.nanoseconds);
    fprintf(stderr, "DEBUG coincidence_window_samples   %d\n", asc->coincidence_window);
    fprintf(stderr, "DEBUG coincidence_sum_threshold_pe %f\n", processor->windowed_sum_threshold_pe);
    fprintf(stderr, "DEBUG sum_threshold_pe             %f\n", processor->sum_threshold_pe);
    fprintf(stderr, "DEBUG enable_muon_coincidence      %d\n", processor->muon_coincidence);

    for (int i = 0; i < asc->ntraces; i++)
      fprintf(stderr,"DEBUG channel %d gain %f threshold %f shaping %d lowpass %f dsp_pre_samples %d dsp_post_samples %d\n",
      asc->tracemap[i], asc->gains[i], asc->thresholds[i], asc->shaping_widths[i], asc->lowpass[i], asc->dsp_pre_samples[i], asc->dsp_post_samples[i]);

  }
  return 1;
}

static inline void event_flag_2char(char* string, size_t strlen, unsigned int event_flags)
{
  assert(strlen >= 9);

  for (size_t i = 0; i < strlen; i++)
    string[i] = '_';

  if (event_flags & EVT_AUX_PULSER)
    string[0] = 'P';
  if (event_flags & EVT_AUX_BASELINE)
    string[1] = 'B';
  if (event_flags & EVT_AUX_MUON)
    string[2] = 'M';
  if (event_flags & EVT_RETRIGGER)
    string[3] = 'R';
  if (event_flags & EVT_EXTENDED)
    string[3] = 'E';

  if (event_flags & EVT_FPGA_MULTIPLICITY)
    string[4] = 'M';
  if (event_flags & EVT_FPGA_MULTIPLICITY_ENERGY_BELOW)
    string[5] = 'L';
  if (event_flags & EVT_FORCE_POST_WINDOW)
    string[6] = '<';
  if (event_flags & EVT_ASUM_MIN_NPE)
    string[7] = '-';
  if (event_flags & EVT_FORCE_PRE_WINDOW)
    string[8] = '>';
  // string[4] = '\0';
}


static inline void st_flag_2char(char* string, size_t strlen, unsigned int st_flags)
{
  assert(strlen >= 5);



  if (st_flags & ST_TRIGGER_FORCE)
    string[0] = 'F';
  if (st_flags & ST_TRIGGER_SIPM_NPE_IN_WINDOW)
    string[1] = 'C';
  if (st_flags & ST_TRIGGER_SIPM_NPE)
    string[2] = 'N';
  if (st_flags & ST_TRIGGER_SIPM_PRESCALED)
    string[3] = 'S';
  if (st_flags & ST_TRIGGER_GE_PRESCALED)
    string[4] = 'G';

  // string[10] = '\0';
}

void LPPFlags2char(LPPState* lpp_state, size_t strlen, char* cstring)
{
  assert(strlen >= 18);

  for (size_t i = 0; i < strlen; i++)
    cstring[i] = '_';

  cstring[0] = lpp_state->write?'W':'D';

  switch(lpp_state->stream_tag) {
    case FCIOConfig: cstring[1] = 'C'; break;
    case FCIOStatus: cstring[1] = 'S'; break;
    case FCIOEvent:  cstring[1] = 'E'; break;
    case FCIOSparseEvent: cstring[1] = 'Z'; break;
    case FCIORecEvent: cstring[1] = 'R'; break;
    default: cstring[1] = '?'; break;
  }
  cstring[2] = '.';

  st_flag_2char(&cstring[3], 5, lpp_state->flags.trigger);
  cstring[8] = '.';

  event_flag_2char(&cstring[9], 9, lpp_state->flags.event);

}
