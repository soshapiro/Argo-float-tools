%put the fill argos in a similar structure as the old struct
%cor files are the ones with GL_ files instead of subfolders and profiles
%preferred method !
clear 

fnams = dir('~/Documents/Research/Chapter2/Argo/NewArgo/DataSelection_NEWALL/*.nc');

pr_grid = [0:2:2000];
sg1_grid = [28:0.01:32.6];
gamma_grid = [24:0.01:28.3];

var=[{'ts'}, {'sp'},{'sa'}, {'tc'}, {'pr'}, {'sg1'}, {'gamma'}];
int = [{'pr'}, {'sg1'}]; %do we want to interpolate on neutral density as well??
const = [{'lat'}, {'lon'}, {'juldDPS'}, {'WMOid'}]; %these variables will be copied over uninterpolated
oxlog = true(length(fnams), 2);
%%

for i = 1:length(fnams)
    i
    c = fnams(i).name;
    
    info = ncinfo(c);
    atts = info.Attributes;
    % find the data mode in the info file
    for j = length(atts):-1:1
        if contains(atts(j).Name, 'mode')
            modej = j;
            break %once you've found it, skip the rest of js
        else
            modej = NaN;
        end
    end %only use D (delayed mode) data
    if isfinite(modej)
        if contains(atts(modej).Value, 'R', 'IgnoreCase', true)
            continue %skip data files that are in realtime mode
        end
    end
    
    %latitude and longitude, important vars
    D(i).WMOid = c(end-9:end-3);
    
    time = ncread(c, 'TIME');
    D(i).juld_qc = ncread(c,'TIME_QC'); %max(D(i).juld_qc) %if you're getting 4s here you should do something
    time(D(i).juld_qc~=1)=nan;
    D(i).juldDPS = days(time)+datetime(1950, 1,1, 0, 0, 0);
    
    D(i).lon = ncread(c, 'LONGITUDE');
    D(i).lat = ncread(c, 'LATITUDE');
    
    D(i).pos_qc=ncread(c, 'POSITION_QC');%get rid of 4;
    %D(i).pos_qc=str2num(ncread(c, 'POSITION_QC'));%get rid of 4
    D(i).lon(D(i).pos_qc~=1)=nan;
    D(i).lat(D(i).pos_qc~=1)=nan;
    
    %add pressure, filter out a gc of 4
    D(i).pr = ncread(c, 'PRES_ADJUSTED');
    pr_qc=ncread(c,'PRES_ADJUSTED_QC'); %1= good 4=bad
    for j=1:size(pr_qc, 2) %for each column of data
        ii=strfind(pr_qc(:, j)',' ');pr_qc(ii,j)='0'; %set blanks (at depth) to 0
        %D(i).pr_qc(:, j)=str2num(pr_qc(:, j)); %convert string to columned nums
        D(i).pr_qc(:, j)=(pr_qc(:, j)); %convert string to columned nums

    end
    D(i).pr(D(i).pr_qc~=1)= nan; %set values with qc of 4 to nan
    D(i).pr(D(i).pr<0)=nan; %set values with negative pr to nan

    %add temp, filter out qc of 4
    try
        D(i).ts = ncread(c, 'TEMP_ADJUSTED');
    catch 
        D(i) = [];
        continue
    end
    ts_qc=ncread(c,'TEMP_ADJUSTED_QC'); %1= good 4=bad
    for j=1:size(ts_qc, 2) %for each column of data
        ii=strfind(ts_qc(:, j)',' ');ts_qc(ii,j)='0'; %set blanks (at depth) to 0
        %D(i).ts_qc(:, j)=str2num(ts_qc(:, j)); %convert string to columned nums
        D(i).ts_qc(:, j)=(ts_qc(:, j)); %convert string to columned nums
    end
    D(i).ts(D(i).ts_qc~=1)= nan; %set values with qc of 4 to nan

    %add salinity, filter out qc of 4
    try
        D(i).sp = ncread(c, 'PSAL_ADJUSTED');        
    sp_qc=ncread(c,'PSAL_ADJUSTED_QC'); %1= good 4=bad
    for j=1:size(sp_qc, 2) %for each column of data
        ii=strfind(sp_qc(:, j)',' ');sp_qc(ii,j)='0'; %set blanks (at depth) to 0
        %D(i).sp_qc(:, j)=str2num(sa_qc(:, j)); %convert string to columned nums
        D(i).sp_qc(:, j)=(sp_qc(:, j)); %convert string to columned nums    
    end
    D(i).sp(D(i).sp_qc~=1)= nan; %set values with qc of 4 to nan
    D(i).sp(abs(D(i).sp-0)<0.05) = nan; %set zero values to nan (some use 0 as a fill value)
    catch 
        D(i).sp = nan(size(D(i).ts));
    end
        %loop through each of the variables I'm interested in
    variables = {'DOXY_ADJUSTED', 'DOX2_ADJUSTED'};
    newvars = {'ox', 'ox2'};
    
    for v = 1:2
        %grab the relevant data from the float
        try %do this if the variable was measured
            str = [newvars{v} '= ncread(c, "' variables{v} '");'];
            eval(str);
            %filter out data with bad QC flags
            str = [variables{v} '_QC']; 
            qc=ncread(c, str); %1= good 4=bad
            for j=1:size(qc, 2) %for each column of data
                ii=strfind(qc(:, j)',' ');qc(ii,j)='0'; %set blanks (at depth) to 0
                str = [newvars{v} '_qc(:, j)=qc(:, j);']; %convert string to columned nums
                eval(str);
            end
            str = [newvars{v} '(' newvars{v} '_qc==4)= nan;']; %set values with qc of 4 to nan
            eval(str);
            %if all of the values are nan, go to catch

        catch %do this if the variable wasn't measured
            str = [newvars{v} '= [NaN(size(D(i).pr))];'];
            eval(str);
            oxlog(i, v) = false;
        end
    end
            


    %calculate the conductivity from gsw
    %adjust the conductivity  . . . 
    %re-compute practical salinity
    
    D(i).sa = gsw_SA_from_SP(D(i).sp,D(i).pr,D(i).lon,D(i).lat); %calculate Absolute SA from practical SA
    D(i).tc = gsw_CT_from_t(D(i).sa,D(i).ts,D(i).pr); %Conservative temp from in-situ temp
    D(i).sg1 = gsw_sigma1(D(i).sa, D(i).tc); %potential density @1000 anomaly from absolute SA and CT
    lat = repmat(D(i).lat', size(D(i).sa, 1), 1);
    lon = repmat(D(i).lon', size(D(i).sa, 1), 1);
    [D(i).gamma, ~, ~] = eos80_legacy_gamma_n(D(i).sp,D(i).ts,D(i).pr, lon, lat);

    % grid the files
    for g = 1:length(int)   
        for c=1:length(const)    %assign the constant vars to each dataset 
            %assign the variable to the new structure eg D_pr.lat = D_ctd(1).lat
            str = ['D_' int{g} '(i).' const{c} '= D(i).' const{c} ';'];
            eval(str);
        end 
        for v = 1:length(var)  %interpolate each of the vars onto each grid
            dat = eval(['D(i).' var{v}]);
            grid = eval([int{g} '_grid']);
            locvar = eval(['D(i).' int{g}]);
            eval(['D_' int{g} '(i).' var{v} ' = profinterp(dat, locvar, grid, "pchip");']);
        end
    end
   

    %calculate PV and add it to the data files
    D_pr(i).PV = PVf(D_pr(i).pr, D_pr(i).sa, D_pr(i).tc, D_pr(i).lat, pr_grid, 41);
    %interp PV onto sg1 grid
    D_sg1(i).PV = profinterp(D_pr(i).PV, D_pr(i).sg1, sg1_grid, "pchip");
    
    %interpolate onto neutral density surfaces
    [D_gamma(i).sp,D_gamma(i).ts,D_gamma(i).pr] = eos80_legacy_neutral_surfaces(D(i).sp,D(i).ts,D(i).pr,D(i).gamma,gamma_grid);
     D_gamma(i).sa = gsw_SA_from_SP(D_gamma(i).sp,D_gamma(i).pr,D(i).lon,D(i).lat);
     D_gamma(i).tc = gsw_CT_from_t(D_gamma(i).sa,D_gamma(i).ts,D_gamma(i).pr);
     D_gamma(i).th = gsw_pt0_from_t(D_gamma(i).sa,D_gamma(i).ts,D_gamma(i).pr); %calculate ptemp from temp, Absolute SA, reference level = 0
     D_gamma(i).PV = profinterp(D_pr(i).PV, D_pr(i).gamma, gamma_grid, "pchip");
end
%remove any empty columns from the structure 
%(if we skipped that i, get rid of it)
D(arrayfun(@(S) isempty(S.lon), D)) = [];
D_pr(arrayfun(@(S) isempty(S.lon), D_pr)) = [];
D_gamma(arrayfun(@(S) isempty(S.ts), D_gamma)) = [];
D_sg1(arrayfun(@(S) isempty(S.lon), D_sg1)) = [];

%%
%{.
lat = vertcat(D_pr(:).lat);
lon = vertcat(D_pr(:).lon);
date = vertcat(D_pr(:).juldDPS);
maxpr = arrayfun(@(S) max(S.pr), D_pr, 'UniformOutput', false);
maxpr = horzcat(maxpr{:})';
maxsg1 = arrayfun(@(S) max(S.sg1), D_sg1, 'UniformOutput', false);
maxsg1 = horzcat(maxsg1{:})';
maxg = arrayfun(@(S) max(S.gamma), D, 'UniformOutput', false);
maxg = horzcat(maxg{:})';
%type1 = arrayfun(@(S) repmat(S.WMOid', 1, length(S.lat))', D_pr, 'UniformOutput', false);
%type = num2cell(vertcat(type1{:}));
type1 = arrayfun(@(S) repmat(str2num(S.WMOid)', 1, length(S.lat)), D_pr, 'UniformOutput', false);
type = num2cell(horzcat(type1{:}))';


properties = [{'pr'}, {'ts'}, {'sa'}, {'tc'}, {'PV'}, {'sg1'}, {'gamma'}];
Pmat = nan(length(pr_grid), length(lat), length(properties));
sg1mat = nan(length(sg1_grid), length(lat), length(properties));
gmat = nan(length(gamma_grid), length(lat), 5);

for p = 1:length(properties)
    Pmat(:, :, p) = horzcat(D_pr(:).(properties{p}));
    sg1mat(:, :, p) = horzcat(D_sg1(:).(properties{p}));
end
p = [1 2 3 4 5];
for i = 1:5
    gmat(:, :, p(i)) = horzcat(D_gamma(:).(properties{p(i)}));
end
%}
%Pmat(:, :, 8) = gsw_geo_strf_dyn_height(Pmat(:, :, 3),Pmat(:, :, 4),Pmat(:, :, 1),p_ref);
%properties(8) = {'dyn_height'};


%Pmat(:, :, 8) = nan(size(Pmat(:, :, 5)));
%sg1mat(:, :, 8) = nan(size(sg1mat(:, :, 5)));
oxlog1 = or(oxlog(:, 1), oxlog(:, 2));

%%
save '~/Documents/Research/Chapter2/hotspotbig.mat' D D_pr D_sg1 D_gamma pr_grid sg1_grid gamma_grid oxlog1
%%
%{.
keep properties Pmat sg1mat gmat lat lon date type maxpr maxsg1 maxg pr_grid sg1_grid gamma_grid
save '~/Documents/Research/Chapter2/Argo/hotspotbigP.mat'
%}
%% 
%{
%add the other core data
%C = load("~/Documents/Research/Data_Files/depthgriddeddatasmall.mat");
C = load("~/Documents/Research/Chapter2/Argo/hotspotbigP.mat");

names = fieldnames(C); 
%concatenate new grids onto old grids
for n = 1:length(names)
    if contains(names{n}, 'mat')
        str = [names{n} '= horzcat(', names{n} ', C.' names{n} ');'];
        eval(str);
    elseif ~contains(names{n}, 'grid') & ~contains(names{n}, 'properties')
        str = [names{n} '= vertcat(' names{n} ', C.' names{n} ');'];
        eval(str);
    end
end
save '~/Documents/Research/Data_Files/hotspotbigP.mat'

%}

%%

