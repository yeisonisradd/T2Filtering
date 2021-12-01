%% READLISTDATA_YW     Reads a Philips .LIST/.DATA file pair
%
% [DATA,INFO] = READLISTDATA_YW(FILENAME)
%
%   FILENAME is a string containing a file prefix or name of the .LIST
%   header file or .DATA binary file, e.g. RAW_001 or RAW_001.LIST or RAW_001.DATA
%
%   DATA is an N-dimensional array holding the raw binary data.
%
%   INFO is a structure containing details from the .LIST header file
%
% Example:
%  [data,info] = readListData('raw_777.list');
%
%  See also: LOADLABRAW
%
%  Dependencies: none
%

%% Revision History
% * 2006.01.01  initial version - brianwelch
% * 2012.08.20  Ananth Madhuranthakam - added extr1 for ASL and ky output 
% * 2018.08.30  Yiming Wang - added T2 correction filter
% * 2021.10.19 Yeison Rodriguez - Updated program to ask user for desired
% T2 filter cut off and to call program for T2 correction filtering
% Program updated to filter the average image instead of individual
% dynamics

function [data,info,ky_range] = readListData_YRV3(name)
%hard code name for testing
%YW 08-31-2018 - start
T2_filter = 1;
%YW 08-31-2018 - end


%obtain the prefix from the file name entered
toks = regexp(name,'^(.*?)(\.list|\.data)?$','tokens');
prefix = toks{1}{1};

%create the name of the files (list or data)
listname = sprintf('%s.list',prefix);
dataname = sprintf('%s.data',prefix);

%open the files for reading
fid = fopen(listname,'r');
if fid~=-1,
    listtext = fread(fid,inf,'uint8=>char')';
    fclose(fid);
else
    error( sprintf('cannot open %s for reading', listname) );
end
%list of the atributes of the complex vector types
attr = {'mix','dyn','card','echo','loca','chan','extr1','extr2','ky','kz','n.a.','aver','sign','rf','grad','enc','rtop','rr','size','offset'};

% clean attr names (will be used as variablenames and fieldnames)
for k=1:length(attr),
    attr{k} = cleanFieldname( attr{k} );
end

pattern = '(?<typ>\w+)';
for k=1:length(attr),
    pattern = sprintf('%s\\s+(?<%s>-?\\d+)',pattern,attr{k});
end
%create a structured cell array of each name token in each match
info = regexp(listtext, pattern, 'names');
%strcmp will compare the two strings and their equivalency
%find will find the index at which this equivalency is true
idxSTD = find( strcmp({info(:).typ},'STD') );

for k=1:length(attr),
    eval( sprintf('%s = sort( str2num( char( unique( {info(idxSTD).%s} ) ) ) );',attr{k},attr{k}) );
end

% DATA is a multi-dimensional array organized as
order = {'kx','ky','kz','loca','dyn','card','echo','mix','aver','extr1','chan'};
%create temporary string for storage
tmpstr = '';

for k=1:length(order),
    %if the string is kx
    
    if strcmp(order{k},'kx')==1,
        tmpstr = sprintf('%s max(size)/4/2',tmpstr);
    else
        %tmpstr = sprintf('%s length(%s)',tmpstr,order{k});
        tmpstr = sprintf('%s max(%s)-min(%s)+1',tmpstr,order{k},order{k});
    end
end
%Data variable is defined
eval( sprintf('data = zeros([%s]);', tmpstr) );

% Ananth
ky_range = zeros(1,2);
ky_range(1) = min(ky);
ky_range(2) = max(ky);

% kx_range = zeros(1,2);
% kx_range(1) = min(kx);
% kx_range(2) = max(kx);

for k=1:length(order),
    %the string is not kx
    if strcmp(order{k},'kx')==0,
        eval( sprintf('tmp = 1 - min(%s(:));',order{k}) );
        eval( sprintf('%s(:) = %s(:) + tmp;',order{k},order{k}) );
        for j=1:length(idxSTD),
            info(idxSTD(j)).(order{k}) = str2double( info(idxSTD(j)).(order{k}) ) + tmp;
        end
    end
end



fid = fopen(dataname,'r','ieee-le');
if fid==-1,
    error( sprintf('cannot open %s for reading', listname) );
end

hwait = waitbar(0,'=========================================================================================');
set( get( findobj(hwait,'type','axes'),'Title') ,'Interpreter','none');
set( get( findobj(hwait,'type','axes'),'Title') ,'String',sprintf('Reading raw data from %s ...', dataname) );

N = length(idxSTD);

%YW 08-15-2018
Echo_space = 3; %supposedly 
T2_GM = 99;   %Stanisz-2005-MRM
T2_Blood = 275; %T2 for kidney's was said to be 275
ch_num = length(chan);
rf_num = length(rf);
%end YW


%%% YR 11.30.2021 - start
size_data = size(data);
length_datasize = length(size_data);

nx_oversample = 2;
nx      = size_data(1)/nx_oversample;
if (length_datasize > 1)
    ny      = size_data(2);
else
    ny = 1;
end
if (length_datasize > 2)
    nz      = size_data(3);
else
    nz = 1;
end
if (length_datasize > 3)
    nlocs   = size_data(4);
else
    nlocs = 1;
end
if (length_datasize > 4)
    ndynamics  = size_data(5);
else
    ndynamics = 1;
end
if (length_datasize > 5)
    ncards  = size_data(6);
else
    ncards = 1;
end
if (length_datasize > 6)
    nechoes = size_data(7);
else
    nechoes = 1;
end
if (length_datasize > 7)
    nmix    = size_data(8);
else
    nmix = 1;
end
if (length_datasize > 8)
    navgs   = size_data(9);
else
    navgs = 1;
end
if (length_datasize > 9)
    lbl_ctrl= size_data(10);
else
    lbl_ctrl = 1;
end
if (length_datasize > 10)
    ncoils  = size_data(11);
else
    ncoils = 1;
end




% Average control data and label data for each coil
ctrl_data = zeros(nx*nx_oversample,ny,nz,ncoils);
lbl_data = zeros(nx*nx_oversample,ny,nz,ncoils);

    tmpDataC = (mean(data(:,:,:,:,:,:,:,:,:,1,:),5));
    tmpDataL = (mean(data(:,:,:,:,:,:,:,:,:,2,:),5));
    for j = 1:ndynamics
         data(:,:,1,1,j,1,1,1,1,1,:) = tmpDataC;
         data(:,:,1,1,j,1,1,1,1,2,:) = tmpDataL;
    end
    size_Test = size(data);

    
    
       
%%% YR 10.19.2021 - start
if T2_filter == 1
    %Call filter program
    [T2_cf_tp] = T2_Filter(mod(ceil(N/ch_num) - 1,rf_num)+1);
    %Save magnitude of data
    T2_cf = abs(T2_cf_tp);
end
%YR 10.19.2021 - end

% ETL = 120 for 3D(Not 140?), 83 for raw005, 81 for raw007, and 81 for raw014
for n=1:N,
    if( fseek(fid, str2num( info(idxSTD(n)).offset ) ,'bof') == 0),
        tmpdata = fread(fid, str2num( info(idxSTD(n)).size )/4 ,'float32');
        %turns into complex double
        tmpdata = tmpdata(1:2:end) + i*tmpdata(2:2:end);
        tmpdata = tmpdata * str2num( info(idxSTD(n)).sign );
        
        %YR
        if T2_filter == 1   
            
            echo_idx = mod(ceil(n/ch_num) - 1,rf_num);
            %It seems here is where the Filter is applied.
            %Note this is data from info not data var.
            tmpdata = tmpdata*T2_cf(echo_idx+1);
        end
        %end YW
        tmpstr='';
        for k=1:length(order),
                %basically only for the first step of the data
                if strcmp(order{k},'kx')==1,
                    tmpstr = sprintf('%s,1:%d', tmpstr, length(tmpdata));
                %all other data goes through here
                else
                    %eval( sprintf('idx = find( %s==str2num(info(idxSTD(n)).%s) );', order{k}, order{k} ) );
                    idx = info(idxSTD(n)).(order{k});
                    tmpstr = sprintf('%s,%d', tmpstr, idx);
                end        
        end
                
        tmpstr(1)=[]; % Delete initial comma'
        
        %tmp data from info var is written to data, passed to main program
        eval( sprintf('data(%s) = tmpdata;', tmpstr) );
        
    else
        error('Cannot FSEEK to offset=%d in data file %s', info(idxSTD(k)).offset,dataname); 
    end
    if mod(n,100)==99,
        waitbar(n/N,hwait);
    end
end

%%
% Squeeze data and the end result is:
% 2D: kx-ky-ndynamics-lbl_ctrl-ncoils (Need to modify for >1 slice)
% 3D: kx-ky-kz-lbl_ctrl-ncoils

fclose(fid);
close(hwait);

% # Complex data vector types:
% # --------------------------
% # STD = Standard data vector (image data or spectroscopy data)
% # REJ = Rejected standard data vector
% #       (only for scans with arrhythmia rejection)
% # PHX = Correction data vector for EPI/GraSE phase correction
% # FRX = Correction data vector for frequency spectrum correction
% # NOI = Preparation phase data vector for noise determination
% # NAV = Phase navigator data vector
% #
% # Other attributes of complex data vectors:
% # -----------------------------------------
% # mix    = mixed sequence number
% # dyn    = dynamic scan number
% # card   = cardiac phase number
% # echo   = echo number
% # loca   = location number
% # chan   = synco channel number
% # extr1  = extra attribute 1 (semantics depend on type of scan)
% # extr2  = extra attribute 2 (   ''       ''   ''  ''  ''  '' )
% # kx,ky  = k-space coordinates in 1st and 2nd preparation direction (spectroscopy data)
% # ky,kz  = k-space coordinates in 1st and 2nd preparation direction (image data)
% # aver   = sequence number of this signal average
% # sign   = sign of measurement gradient used for this data vector (1 = positive, -1 = negative)
% # rf     = sequence number of this rf echo (only for TSE, TFE, GraSE)
% # grad   = sequence number of this gradient echo (only for EPI/GraSE)
% # enc    = encoding time (only for EPI/GraSE)
% # rtop   = R-top offset in ms
% # rr     = RR-interval length in ms
% # size   = data vector size   in number of bytes (1 complex element = 2 floats = 8 bytes)
% # offset = data vector offset in number of bytes (first data vector starts at offset 0)
% #
% # The complex data vectors are represented as binary data in little endian single precision IEEE float format.
% #
% # Please note that complex data vector attributes which are not relevant for a certain type of vector
% # may have arbitrary values!

% # Identifying attributes of complex data vectors:
% # -----------------------------------------------
% # The next table specifies the identifying attributes for each type of complex data vector:
% #
% # typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    aver  sign  rf    grad  enc   rtop  rr    size   offset
% # --- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ------ ------
% #
% # STD   *     *     *     *     *     *     *     *     *     *     *     *     *     *
% # REJ   *     *     *     *     *     *     *     *     *     *     *     *     *     *
% # PHX   *                 *     *     *                                   *     *     *
% # FRX   *                 *     *     *                                   *            
% # NOI                           *     *                                                
% # NAV   *     *     *     *     *     *     *     *     *     *     *     *     *     *

