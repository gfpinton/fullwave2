function [apmlx1 apmlu1] = pml_gradient_mask2(nX,nY,nbdy,M,crd,apmlx1,apmlu1)


pmlmask=zeros(nX,nY);
if(crd==1)
    for i=1:nbdy
        pmlmask(i+(nX-M-nbdy+1),:)=(i/nbdy);
        pmlmask(M+nbdy+1-i,:)=(i/nbdy);
    end
    pmlmask(1:M,:)=1; pmlmask(nX-M+1:nX,:)=1;
end

if(crd==2)
    for i=1:nbdy
        pmlmask(:,i+(nY-M-nbdy+1))=(i/nbdy);
        pmlmask(:,M+nbdy+1-i)=(i/nbdy);
    end
    pmlmask(:,1:M)=1; pmlmask(:,nY-M+1:nY,:)=1;
end


tmpx=apmlx1;
tmpu=apmlu1;

apmlx1=(tmpx.*(1-pmlmask)+tmpu.*pmlmask);
apmlu1=(tmpu.*(1-pmlmask)+tmpx.*pmlmask);

