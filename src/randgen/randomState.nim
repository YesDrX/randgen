import locks
import threadpool

when defined(js):
  type RandU = uint32
  const randMax = 4_294_967_295u32
else:
  type RandU = uint64
  const randMax = 18_446_744_073_709_551_615u64

import random
import tables

when defined(js):
  type RandInt = uint32
else:
  type RandInt = uint64

type
  RandState* = object
    a0, a1 : RandInt
  
var
  states : Table[int, Rand]
  thread_ids : seq[int]
  lock : Lock
  

proc getThreadIDforRandState(): int =
  when compileOption("threads"):
    result = getThreadID()
  else:
    result = 0

proc initRandState*(seed: SomeInteger) =
  let
    threadId = getThreadIDforRandState()
  {.gcsafe.}:
    acquire(lock)
    if not thread_ids.contains(threadId):
      thread_ids.add(threadId)
      states[threadId] = initRand(seed xor thread_ids.len)
    release(lock)

proc nextState*(): uint64 =
  let
    threadId = getThreadIDforRandState()
  {.gcsafe.}:
    result = states[threadId].next()

proc work() =
  {.gcsafe.}:
    initRandState(100)
    echo states
    discard nextState()
    echo states

when isMainModule:
  for i in 0 .. 5:
    spawn work()
  sync()