function [im, SupNeighborTable] = merge_superpixels(im, param)

%This converts a non_ordered and non_connected superpixel image into 
% an ordered superpixel image
%
% inputs:
% im : N by M matrix depends on the image size       
% SmallThre : the samllest Sup size
%
% outputs:
% im_order : N by M matrix which is ordered

%%%
% For each of the superpixel indicies, find all of the
% disconnected components with that label.  if the component is
% small (< 200 pixels) or isn't the biggest component with that
% index, call analysesupinpatch with the outline of the component
% to replace it's label with the one that is most common in the outline
%%%

if nargin <2
   param.SmallThre = 5; % smallest sup size
end
SupNeighborTableFlag = 1;

[yn, xn] = size(im);
NuSup = unique(im(:))';
SE = strel('octagon',3);  
for i=NuSup

        % label connected component
        temp = zeros(size(im));
        temp(im(:,:)==i)=1;
        [L,num] = bwlabel(temp,4);

        % find the main piece
         [maxL, dum]= mode(L(L~=0));
%        his = histc(L(:), 1:num);
%        [dum maxL ]= max(his);

           if dum > param.SmallThre;
              SupMerg = setdiff(1:num,maxL);
           else
              SupMerg = 1:num;
           end

           for k = SupMerg
               mask = L==k;
               % then assign those pixels to mostlikely 3 by 3 neighborhood
                mask_dilate = imdilate(mask,SE);
%                mask_dilate = mask | [zeros(yn,1) mask(:,1:(end-1))] ...
%                                   | [mask(:,2:(end)) zeros(yn,1)] ...
%                                   | [zeros(1,xn) ;mask(1:(end-1),:)] ...
%                                   | [mask(2:(end),:); zeros(1, xn)] ....
%                                   | [[zeros(yn-1,1) mask(2:end,1:(end-1))]; zeros(1,xn)]...
%                                   | [[mask(2:end,2:(end)) zeros(yn-1,1) ]; zeros(1,xn)]...
%                                   | [zeros(1,xn); [zeros(yn-1,1) mask(1:(end-1),1:(end-1))]]...
%                                   | [zeros(1,xn); [mask(1:(end-1),2:(end)) zeros(yn-1,1)]]...
%                                   ;                          
               mask_dilate(mask) = 0;
%                im(mask) = analysesupinpatch(im(mask_dilate));%hard work
               im(mask) = mode(im(mask_dilate));
           end
end

% merge the small superpixel with the surrrounding one if it's neighbor is only one
 MaxSupIndex = max(NuSup(:));
 SupNeighborTable = sparse(MaxSupIndex,MaxSupIndex);

if SupNeighborTableFlag
   for i = 1:((xn-1)*yn)
     
       % registed the neoghbor in right
       SupNeighborTable(im(i),im(i+yn)) = 1;
       SupNeighborTable(im(i+yn),im(i)) = 1;

       % registed the neoghbor in below
       if mod(i,yn) ~=0
          SupNeighborTable(im(i),im(i+1)) = 1;
          SupNeighborTable(im(i+1),im(i)) = 1;
       end
   end

   % find out the single neighbor ones and merge them with neighbors
   SingleTar = sum( SupNeighborTable,1);
   for i = find(SingleTar == 1)
       mask = im == i;
       im(mask) = find(SupNeighborTable(:,i) == 1);
       SupNeighborTable(:,i) = 0;
       SupNeighborTable(i,:) = 0;
   end

end
return;
