import PEATSAmenu
import sys

if len(sys.argv) < 2 or len(sys.argv) > 3:
    print("Syntax:")
    print("python run_locally.py script.emtgopt [script_result.emtg]")
else:
    
    sys.path.append("/Utilities/emtg/PyEMTG/")
    import MissionOptions
    
    MO = MissionOptions.MissionOptions(sys.argv[1])
    
    MO.run_inner_loop = 0
    
    MO.override_working_directory = 0
    
    MO.override_mission_subfolder = 0
    
    MO.quiet_NLP = 0
    
    MO.HardwarePath = "/Utilities/emtg/HardwareModels/"
    
    if len(sys.argv) == 3:
        import Mission
        M = Mission.Mission(sys.argv[2])
        new_dv = []
        for idx in range(len(M.DecisionVector)):
            new_dv.append([M.Xdescriptions[idx],M.DecisionVector[idx]])
        MO.trialX = new_dv
    
    MO.write_options_file(sys.argv[1])
    
    if MO.SpacecraftModelInput == 1:
        import PyHardware
    
        SO = PyHardware.SpacecraftOptions(MO.HardwarePath + "/" + MO.SpacecraftOptionsFile)
    
        for idx in range(SO.getNumberOfStages()):
            StageO = SO.getStageOptions(idx) 
            PropO = StageO.getElectricPropulsionSystemOptions()
            # print PropO.getThrottleTableFile()
            
            # print PropO.getThrottleTableFile()
            StageO.setElectricPropulsionSystemOptions(PropO)
            SO.setStageOptions(idx,StageO)
    
        SO.write_output_file(MO.HardwarePath + "/" + MO.SpacecraftOptionsFile)