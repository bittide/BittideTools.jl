module Framing

export write_popped_data_to_processor_memory, step_processor_one_cycle, read_data_from_processor_memory

# Functions called by the Frames module of Callisto
# and implemented by Elara.
write_popped_data_to_processor_memory(fdp::Nothing, i, port, localtick, data) = return
step_processor_one_cycle(fdp::Nothing, i, localtick) = return
read_data_from_processor_memory(fdp::Nothing, i, port, localtick)::UInt64 = 0x0



end

