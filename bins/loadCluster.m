function [] = loadCluster(nCore,nMem,nQueue,nWall)
%load the cluster environment for corbra jobs
%any input parameter should be in str datatype; e.g., '1024' for memory.
%but core number should be a integer
delete(gcp('nocreate'));%close all existing pools to avoided error
%import profile with a ramdom name to avoid conflict,
%tmp = parallel.importProfile('/home/xl95w/ghpcc_xl_profile.settings')
c = parcluster('ghpcc local R2019a');
c.AdditionalProperties.MemUsage = nMem;
c.AdditionalProperties.WallTime = nWall;
c.AdditionalProperties.QueueName = nQueue;
c.AdditionalProperties.EmailAddress = 'xuhang.li@umassmed.edu';
c.saveProfile
p = parpool(tmp,nCore)
%pctRunOnAll('initCobraToolbox(false)');
%parallel.internal.ui.MatlabProfileManager.removeProfile(tmp);
end
