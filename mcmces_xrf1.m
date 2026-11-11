function [Mchain,logP,rate,initprop]=mcmces_xrf1(Prior,Likelihood,N,Nwalkers,Nkeep,StepSize,ThinChain)
% 
% INPUTS:
%  Prior: prior information (all volume fraction solutions should be <1 & >0
%  Likelihood: the likelihood function
%  N: number of unknowns
%  Nwalkers: numner of walkers (chains, >10)
%  Nstep: number of total recorded steps per walker
%  StepSize: the size of the proposed steps (low maintenance)
%  ThinChain: only record every # steps
%
% OUTPUTS:
%  Mchain: A N*Nwalkers*Nstep matrix Markvo chains (without burn-in)
%  logP: A N*Nwalkers*Nstep matrix of log probabilities for each model in 
%    the models
%  rate: acceptance rate
%  initwalker: initial proposals


% initial walker
initprop=rand(N,Nwalkers);   
initprop = initprop./(sum(initprop,1));

Mchain=nan(N,Nwalkers,Nkeep); 
Mchain(:,:,1)=initprop;

%calculate logP state initial pos of walkers
% tic
logP=zeros(1,Nwalkers,Nkeep);
for wix=1:Nwalkers
    logP(1,wix,1)=Likelihood(initprop(:,wix));   
end

curm=Mchain(:,:,1);
curlogP=logP(:,:,1);
% acceptance
reject=zeros(Nwalkers,1); accept=zeros(Nwalkers,1); total=zeros(Nwalkers,1);

%%
for row=1:Nkeep
%     row
    for jj=1:ThinChain
        
        rix=mod((1:Nwalkers)+floor(rand*(Nwalkers-1)),Nwalkers)+1;
        zz=((StepSize - 1)*rand(1,Nwalkers) + 1).^2/StepSize;
        proposedm=curm(:,rix) - bsxfun(@times,(curm(:,rix)-curm),zz);
        logrand=log(rand(2,Nwalkers));
        %logrand=(rand(2,Nwalkers));
        
        for wix=1:Nwalkers
            
            acceptstep=true;
                                
            if logrand(1,wix) < log(zz(wix)) % zz(wix)%
                             
                if Prior(proposedm(:,wix))
                                        
                    proposedlogP=Likelihood(proposedm(:,wix)); %%
                    
                    if logrand(2,wix) > proposedlogP-curlogP(wix) 
                    %if logrand(2,wix) > proposedlogP / curlogP(wix) 
                        acceptstep=false;                        
                    end
                else
                    acceptstep=false;
                end
                
            else
                acceptstep=false;
            end
            
            if acceptstep
                curm(:,wix)=proposedm(:,wix);
                curlogP(:,wix)=proposedlogP;
                accept(wix)=accept(wix)+1;
            else
                reject(wix)=reject(wix)+1;
            end
            total(wix)=total(wix)+1;
        end

    end
    Mchain(:,:,row)=curm;
    logP(:,:,row)=curlogP;
end
% toc
rate=mean(accept./total);

