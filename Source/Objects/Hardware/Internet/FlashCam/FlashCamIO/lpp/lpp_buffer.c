#include "lpp_buffer.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <lpp_state.h>

LPPBuffer* LPPBufferCreate(unsigned int buffer_depth, Timestamp buffer_window)
{
  LPPBuffer *buffer = (LPPBuffer *) calloc(1, sizeof(LPPBuffer));

  buffer->max_states = buffer_depth + 1;
  buffer->buffer_window = buffer_window;

  buffer->lpp_states = (LPPState*)calloc(buffer->max_states, sizeof(LPPState));

  buffer->insert_state = 0;
  buffer->nrecords_inserted = 0;

  return buffer;
}


void LPPBufferDestroy(LPPBuffer* buffer)
{
  assert(buffer);
  assert(buffer->lpp_states);
  free(buffer->lpp_states);

  free(buffer);

}


LPPState* LPPBufferGetState(LPPBuffer* buffer, int offset)
{
  if (!buffer)
    return NULL;

  int index = (buffer->insert_state + buffer->max_states - 1 + offset) % buffer->max_states;

  if (offset == 0 || offset == 1) {
    return &buffer->lpp_states[index];

  } else if (offset < 0) {
    if (-offset >= buffer->nrecords_inserted || -offset > buffer->max_states - 1) {
      return NULL;
    }
    
    return &buffer->lpp_states[index];

  } else {
    return NULL;

  }
}


LPPState* LPPBufferPeek(LPPBuffer* buffer)
{
  LPPState* lpp_state = LPPBufferGetState(buffer, 1);

  if (lpp_state->in_buffer) {
    /* the event has not been fetched from the buffer */
    return NULL;
  }

  return lpp_state;
}

void LPPBufferCommit(LPPBuffer* buffer)
{
  LPPState* lpp_state = LPPBufferGetState(buffer, 1);

  lpp_state->in_buffer = 1;
  buffer->insert_state = (buffer->insert_state + 1) % buffer->max_states;
  buffer->nrecords_inserted++;
  if (lpp_state->contains_timestamp)
    buffer->buffer_timestamp = timestamp_sub(lpp_state->timestamp,buffer->buffer_window);

  buffer->fill_level++;
}

LPPState* LPPBufferFetch(LPPBuffer* buffer)
{
  LPPState* lpp_state = &buffer->lpp_states[buffer->fetch_state];

  if (lpp_state && lpp_state->in_buffer &&
    (timestamp_greater(buffer->buffer_timestamp, lpp_state->timestamp)|| buffer->flush_buffer)) {

    // advance to the next possible send state
    buffer->fetch_state = (buffer->fetch_state + 1) % buffer->max_states;

    buffer->nrecords_fetched++;
    buffer->fill_level--;

    // record is handed off, forget about it ... until we reuse it.
    lpp_state->in_buffer = 0;

    return lpp_state;
  }
  
  return NULL;
}

int LPPBufferFlush(LPPBuffer* buffer)
{
  buffer->flush_buffer = 1;
  return buffer->fill_level;
}


