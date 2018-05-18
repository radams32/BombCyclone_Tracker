% ========================================================================
%  BOMB CYCLONE TRACKER
%
%  This tracking algorithm uses a 3 step procedure to identify Bombs
%  Cyclones
%       1.) Defines closed low pressure centers
%       2.) Stitches low pressures together as storms through time
%       3.) Identifies from the storms which ones rapidly intensify
%               - this extracts the lat of most rapid intensification
%               - uses pressure deepening metric based off of this lat
%
%  The utility of this tracking algorithm spans any gridded dataset.
%
%  The code was developed to track low pressure in NRR, but can easily be
%  adapted for any gridded reanalysis/GCM
% =========================================================================
%% Closed Lows

lows_ALL = [];
counter = 1;

disp('Finding all low pressure centers...')

for timestamps = 1:1:size(alldata,1);
    %disp(timestamps);
    %Check all grid points for closed low criteria
    for grid_cell = 1:1:size(alldata,2);
        %The only data to loop on are the interior points, leaving 1 grid
        %cell surrounding cells to loop. This checks if it is an
        %interior point
        if ((grid_cell > size(longused,2)+1) && (grid_cell < (size(longused,2)*(size(latused,1)-1))) && (mod(grid_cell, size(longused,2)) > 1));
            %Check all queens case grid points around grid_cell and check
            %if less than or equal to. Closed low
            if ((alldata(timestamps,grid_cell)<= alldata(timestamps, grid_cell-size(longused,2))) && (alldata(timestamps,grid_cell)<= alldata(timestamps, grid_cell+size(longused,2))) && (alldata(timestamps,grid_cell)<= alldata(timestamps, grid_cell-1)) && (alldata(timestamps,grid_cell)<= alldata(timestamps, grid_cell+1)) && (alldata(timestamps,grid_cell)<= alldata(timestamps, grid_cell-size(longused,2)-1)) && (alldata(timestamps,grid_cell)<= alldata(timestamps, grid_cell-size(longused,2)+1)) && (alldata(timestamps,grid_cell)<= alldata(timestamps, grid_cell+size(longused,2)-1)) && (alldata(timestamps,grid_cell)<= alldata(timestamps, grid_cell+size(longused,2)+1)));
                %Document for each closed low
                longitude = longused(1, mod(grid_cell,size(longused,2)));
                latitude = latused(ceil(grid_cell/size(longused,2)),1);
                pressure = alldata(timestamps, grid_cell);
                uni_time = timedata(timestamps,2);
                month = timedata(timestamps,4);
                day = timedata(timestamps,5);
                year = timedata(timestamps,3);
                time = timedata(timestamps,6);
                %Write out first low no matter what
                if counter == 1; 
                    lows_ALL(counter,:) = [latitude, longitude, pressure, uni_time, month, day, year, time];
                    counter = counter + 1;
                else 
                    %Write out additional lows IF they don't share a unique
                    %time OR if they are separated by at least 2.5deg.
                    %This should accept only different low centers.
                    if (uni_time ~= lows_ALL(counter-1,4)) || (abs(latitude - lows_ALL(counter-1,1)) > 2.5) || (abs(longitude - lows_ALL(counter-1,2)) > 2.5);
                        lows_ALL(counter,:) = [latitude, longitude, pressure, uni_time, month, day, year, time];
                        counter = counter + 1;
                    end
                end
                clear longitude latitude pressure uni_time month day year       
            end
        end
    end
end
%% Cyclone Stitching
disp('Stitching low pressures together for storm data...')

cyclone_duration = 1; %this value will count the length of the current low pressure system
DegDis = 7.5; %degree displacement for successive 6hrly timestamps, increments of 2.5deg
cyclone_counter = 1; %number of cyclones. Increments as it finds new cyclones
Cyclones = {}; %cell array that houses all of the cyclones

for lows = 1:1:size(lows_ALL,1)
    cc = lows;
    %disp(lows);
    if cyclone_duration == 1
        Cyclones{cyclone_counter,1}(cyclone_duration,:) = lows_ALL(lows,:);
        cyclone_duration = cyclone_duration + 1;
    end

    if cyclone_duration > 1
        ccc = 1;
        for checks = 1:1:50;
            %disp(checks)
            time1 = lows_ALL(cc,4);
            time2 = lows_ALL(cc+ccc,4);
            %disp([time1, time2])
            lat_disp = lows_ALL(cc+ccc,1)-lows_ALL(cc,1);
            long_disp = lows_ALL(cc+ccc,2)-lows_ALL(cc,2);
            if (time2 - time1 == 0.25) && ((lat_disp <= DegDis && lat_disp >= -2.5) && (long_disp <= DegDis && long_disp >= -2.5));
                Cyclones{cyclone_counter,1}(cyclone_duration,:) = lows_ALL(cc+ccc,:);
                cc = cc+ccc;
                %disp('cc is equal to...')
                %disp(cc)
                cyclone_duration = cyclone_duration + 1;
                ccc = 1;
            else
                ccc = ccc + 1;
            end
        end
        cyclone_counter = cyclone_counter + 1;
        cyclone_duration = 1;
        clear checks
    end
end
%% Extract Bombs
disp('Extracting cyclones that deepen rapidly for their latitude...')

clear cc
press_diffs = [];
Bomb_counter = 1;
Bombs = {};
rapid_int_data = [];

%Cyclone count
for cc = 1:1:size(Cyclones,1)
    %disp('Cyclone count');
    %disp(cc);
    duration = size(Cyclones{cc,1},1);
    %loop for each timestep of each cyclone
    if duration > 1
        for cd = 1:1:duration-1
            %disp('Duration')
            %disp(cd)
            %if the storm lasted less than 4 timesteps then just consider the
            %whole period...else, only consider 24hr windows
            if (duration - cd < 4)
                press_diffs(cd,1) = Cyclones{cc,1}(cd,3) - min(Cyclones{cc,1}(cd+1:end,3));
            else
                press_diffs(cd,1) = Cyclones{cc,1}(cd,3) - min(Cyclones{cc,1}(cd+1:cd+4,3));    
            end
        end
        %find the maximum pressure fall for whole storm and also its index in
        %the variable press_diffs
        [MM,II] = max(press_diffs);
    
        %convert the latitude of the cyclone at the time associated with II
        rad = deg2rad(Cyclones{cc,1}(II,1));
        %disp(press_diffs)
        %disp(press_diffs)
    
        press_diffs = [];
        %calculate if the max pressure fall is greater than Bergeron criteria
        %for explosive cyclogenesis
        if (MM >= (24 * (sin(rad)/sin(sqrt(3)/2))))
            if (Bomb_counter == 1)
                Bombs(Bomb_counter,1) = Cyclones(cc,1);
                rapid_int_data(Bomb_counter,:) = Cyclones{cc,1}(II,:);
                Bomb_counter = Bomb_counter + 1;
            elseif (Bomb_counter > 1)
                dups_check = (Cyclones{cc,1}(1,4) ~= Bombs{Bomb_counter-1,1}(1:end,4));
                if (min(dups_check)~= 0)
                    Bombs(Bomb_counter,1) = Cyclones(cc,1);
                    rapid_int_data(Bomb_counter,:) = Cyclones{cc,1}(II,:);
                    Bomb_counter = Bomb_counter + 1;
                end 
            end
        end
    end
end
