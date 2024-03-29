classdef Animal 
    %ANIMAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ProbeMeta
        Code
        Species
        Strain
        Sex
        Age
        Weight
        GeneticLine
        VirusInjection
        VirusCoordinates
        VirusInjectionDate
        SurgeryDate
        TargetAnatomy
        Anesthesia
        Analgesics
        Antibiotics
        SurgicalComplications
        SurgicalNotes

    end
    
    methods
        function obj = Animal(struct)
            %ANIMAL Construct an instance of this class
            %   Detailed explanation goes here
            obj.Age=struct.Age;
            obj.Analgesics=struct.Analgesics;
            obj.Anesthesia=struct.Analgesics;
            obj.Antibiotics=struct.Antibiotics;
            obj.Code=struct.Code;
            obj.GeneticLine=struct.GeneticLine;
            obj.Sex=struct.Sex;
            obj.Species=struct.Species;
            obj.Strain=struct.Strain;
            obj.SurgeryDate=struct.SurgeryDate;
            obj.SurgicalComplications=struct.SurgicalComplications;
            obj.SurgicalNotes=struct.SurgicalNotes;
            obj.TargetAnatomy=struct.TargetAnatomy;
            obj.VirusCoordinates=struct.VirusCoordinates;
            obj.VirusInjection=struct.VirusInjection;
            obj.VirusInjectionDate=struct.VirusInjectionDate;
            obj.Weight=struct.Weight;
            obj.ProbeMeta=struct.ProbeName;
        end
        function probe=getProbe(obj)
            probename=obj.ProbeMeta;
            sde=experiment.SDExperiment.instance.get;
            probeFolder=sde.FileLocations.General.ProbeFolder;
            list=dir(fullfile(probeFolder,strcat('*',probename,'*')));
            probe=neuro.probe.Probe(fullfile(list.folder,list.name));
        end
    end
end

