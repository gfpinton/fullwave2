function [slicex] = generate_c_scat(c0, cscale, scat_density, scat_size_samp, nY, nZ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% generate a field of random amplitude, random position scatterers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% speed of sound %
%%%%%%%%%%%%%%%%%%
mini_slicex = rand(1,round(nY/scat_size_samp)*round(nZ/scat_size_samp));
zero_vec = zeros(1,length(mini_slicex));
idx = find(mini_slicex>scat_density);
mini_slicex(idx) = zero_vec(idx);
mini_slicex = reshape(mini_slicex,round(nY/scat_size_samp),round(nZ/scat_size_samp))/scat_density;
  
slicex = zeros(nY,nZ);
for j=1:floor(nY/scat_size_samp)
  for k=1:floor(nZ/scat_size_samp)
    for jj=1:scat_size_samp
      for kk=1:scat_size_samp
	slicex((j-1)*scat_size_samp+jj,(k-1)*scat_size_samp+kk) = mini_slicex(j,k);
      end
    end
  end
end
%idz = find(z>scat_depth);
%slicex(:,idz) = zeros(size(slicex(:,idz)));
  
slicex = slicex*c0*cscale;
slicex = slicex+ones(nY,nZ)*c0;
