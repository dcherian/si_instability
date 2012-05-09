function [npvl cpvl cpvr npvr] = find_region(xax,data)
    zmid = ceil(size(data,2)/2);
    xmid = ceil(size(data,1)/2);
    
    ind = find(data(:,zmid) < 0);
     npvl = xax(ind(1)); % negative pv left
     npvr = xax(ind(end)); % and right
     
    ind = find(data(:,zmid) ==  data(xmid,zmid));
     cpvl = xax(ind(1)); % constant pv left
     cpvr = xax(ind(end)); % and right