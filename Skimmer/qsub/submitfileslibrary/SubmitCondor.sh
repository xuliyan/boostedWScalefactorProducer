# Condor submit file generated from template 
# By WTagging scale factors framework
Universe               = vanilla
Getenv                 = True
Executable             = [executable]
#Requirements           = Memory >= 16 GB
Request_memory         = 2 GB
Arguments              = $(arguments) 
Log                    = SingleJobLog$(date)_$(Process).txt
#Input                  = arguments.txt
Output                 = SingleJobSTDOUT$(date)_$(Process).txt
Error                  = SingleJobSTDERR$(date)_$(Process).txt
#Transfer_output_remaps = "output-$(process).txt = results/output-$(process).txt"
transfer_input_files  = $HOME/private/x509up
#file_remaps = "*x509up* = x509up_my"
should_transfer_files = YES
transfer_output_files = $(output) 
when_to_transfer_output = ON_EXIT

+JobFlavour = "[queue]" 
Queue arguments, output from [arguments]  