function [ chanlocs ] = Get_Standard_BESA_Coordinates_From_EEGLAB( strEEGLABToolboxPath )

filename = [strEEGLABToolboxPath,'plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp'];
chanlocs = readlocs(filename);

end