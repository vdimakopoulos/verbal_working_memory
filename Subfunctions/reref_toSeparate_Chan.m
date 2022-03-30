function dataReref = reref_toSeparate_Chan(data, refChan, chansSelection,refmethod)

cfg = [];
cfg.reref = 'yes';
cfg.refmethod = refmethod;
cfg.refchannel = refChan;
cfg.channel  = chansSelection;
dataReref = ft_preprocessing(cfg,data);
end



