function [v,v1,v2,varx,vary,cov] = stochastic_simulation2(N,s,u,steps,...
                                    collect_data,start_time,end_time,outputfile)
% The following source code is a minor modification of code originally made 
% availabe by Pearce MT and Fisher DS, obtained from the site:
% 
% https://datadryad.org/resource/doi:10.5061/dryad.f36v6 
%
% Stochastic simulation for two chromosome model. 
% Each step is one generation. Rates (such as s,u,r) are per generation. 
% The expected size of each subpop after selection, mutation, and mating is computed.
% If expected size < :cutoff: then the size is drawn from Poisson
% distribution with mean equal to the expected size. Otherwise the
% subpopulation size = expected size. 

% output :v, v1, v2:        total rate of adaptation, and rate of adaptation in traits 1 and 2. 
% output :varx, vary, cov:  time averaged variances and covariance in traits 1 and 2.

% input :N: population size
% input :s: effect size of beneficial mutation. 
% input :u: mutation rate per chromosome
% input :steps: number of steps for simulation.
% input :collect_data: true/false - collect detailed data on 2d distr. per generation
% input :start_time: start time for collecting detailed data on 2d distribution
% input :end_time: end time for collecting detailed data on 2d distribution
% input :outputfile: string with filename where detailed data will be stored

% initialization of variables
pop=N;              %abundances of a classes
fit=0;              %total fitness of a class
fitx=0;             %fitness in trait 1 of a class
fity=0;             %fitness in trait 2 of a class
nosefitness = 0;    %total fitness of the front

meanfitness = 0;    %mean fitness of the population
meanfitx = 0;       %mean fitness in x 
meanfity = 0;       %mean fitness in y

varx = 0;           %variance in trait 1
vary = 0;           %variance in trait 2
cov = 0;            %covariance between trait 1 and 2

Na = 0;             %actual population size

cutoff=10/s;        %population cutoff for stochasticity


if (collect_data)   %store parameters used in simulation
    fileID = fopen([outputfile '-0.txt'],'w');
    fprintf(fileID,'%f,%f,%f',N, s, u);
    fclose(fileID);
end

fileID1 = fopen([outputfile '-1.txt'],'w'); %file for all other 2d wave data per generation
fileID2 = fopen([outputfile '-2.txt'],'w'); %file for data on classes per generation
fileID3 = fopen([outputfile '-3.txt'],'w'); %file for data on abundances per generation

for timestep=1:steps   
        
    %%%%%%%%%%%%%%%%%%%
    % Remove columns of zeros
    while any(pop(:,1))==0
        pop(:,1)=[]; 
        fit(:,1)=[];
        fity(1)=[];
    end
    while any(pop(1,:))==0 
        pop(1,:)=[];
        fit(1,:)=[];
        fitx(1)=[];
    end
    
    % Add columns for padding
    dim=size(pop); % stores size of array containing population abundances
    if any(pop(:,dim(2)))==1  % checks for expansion of the array (front) 
        pop(:,dim(2)+1)=zeros(dim(1),1);
        fit(:,dim(2)+1)=fit(:,dim(2))+ones(dim(1),1);
        fity(dim(2)+1)=fity(dim(2))+1;
    end
    
    dim=size(pop);
    if any(pop(dim(1),:))==1
        pop(dim(1)+1,:)=zeros(1,dim(2));
        fit(dim(1)+1,:)=fit(dim(1),:)+ones(1,dim(2));
        fitx(dim(1)+1)=fitx(dim(1))+1;
    end
    
    %%%%%%%%%%%% 
    % Find expected frequencies after selection, mutation, and
    % recombination
    dim=size(pop);
    meanfit=sum(sum(times(pop,fit)))/N;
    
    freq=pop/N;
    newfreq=times(exp(s*(fit-meanfit)),freq); %after selection
    newfreq=newfreq/sum(sum(newfreq)); %make sure frequencies still add to one.
    
    z1=zeros(dim(1),1);
    mutatex=[z1 newfreq];
    mutatex(:,dim(2)+1)=[]; % newfreq already has padding, get rid of extra padding from shift
    z2=zeros(1,dim(2));
    mutatey=[z2; newfreq];
    mutatey(dim(1)+1,:)=[]; % newfreq already has padding, get rid of extra padding from shift
    
    nomutate=(1-2*u)*newfreq;
    postmutate=nomutate+(u)*mutatex+(u)*mutatey;
    
    noeventfreq=nomutate; %no mutation or recombination;
    mutationfreq=((u)*mutatex+(u)*mutatey); %mutation but no recombination
        
    newfreq=noeventfreq+mutationfreq;
    
    % For subpopulations with size less than the stoch_cutoff, draw size
    % from poisson distribution. Otherwise, size = N*(expected frequency).
    
    newpop=N*newfreq;
    stoch=newpop<cutoff;
    stochpop=poissrnd(newpop(stoch)); %does poisson random number for fitness classes below stochasticity cutoff
    newpop(stoch)=stochpop;
    
    newpop=round(newpop);    
    Na = sum(sum(newpop));
    
    meanfitness=sum(sum(times(newpop,fit)))/Na;
    meanfitx = sum(sum(newpop,2).*fitx')/Na;
    meanfity = sum(sum(newpop,1).*fity)/Na;
    
    nosefitness = max(fit(newpop>(1/s)./(fit-meanfitness)));
    
    fitx_arry = fitx'*ones(1,dim(2));
    fity_arry = ones(dim(1),1)*fity;
    
    % recompute time-average of variances and covariances
    if timestep == 1
        varx = sum(sum(times(newpop,(fitx_arry-meanfitx).^2)))/Na;
        vary = sum(sum(times(newpop,(fity_arry-meanfity).^2)))/Na;
        cov = sum(sum(times(newpop,(fitx_arry-meanfitx).*(fity_arry-meanfity(timestep)))))/Na;
    else
        varx = (1/timestep)*((timestep-1)*varx + sum(sum(times(newpop,(fitx_arry-meanfitx).^2)))/Na);
        vary = (1/timestep)*((timestep-1)*vary + sum(sum(times(newpop,(fity_arry-meanfity).^2)))/Na);
        cov = (1/timestep)*((timestep-1)*cov + sum(sum(times(newpop,(fitx_arry-meanfitx).*(fity_arry-meanfity))))/Na);
    end
    
    pop=newpop;
    
    if( collect_data && (timestep >= start_time) && (timestep <= end_time) )
        
        % compute the covariance using only classes at the front
        indx_front = (fit==nosefitness);
        meanfitx_front = sum(pop(indx_front).*fitx_arry(indx_front))/sum(pop(indx_front));
        meanfity_front = sum(pop(indx_front).*fity_arry(indx_front))/sum(pop(indx_front));
        front_cov = sum(pop(indx_front).*(fitx_arry(indx_front)-meanfitx_front).*(fity_arry(indx_front)-meanfity_front))/sum(pop(indx_front));
        
        % compute variances, covarainces and population load
        sigmax2 = sum(sum(times(newpop,(fitx_arry-meanfitx).^2)))/Na;
        sigmay2 = sum(sum(times(newpop,(fity_arry-meanfity).^2)))/Na;
        sigmaxy = sum(sum(times(newpop,(fitx_arry-meanfitx).*(fity_arry-meanfity))))/Na;
        pop_load = nosefitness - meanfitness;
        
        % compute v1, v2, sigma12, G eigenvalues and orientation 
        [D,L] = eig([sigmax2 sigmaxy; sigmaxy sigmay2]);
        evec1 = [cosd(45) sind(45); -sind(45) cosd(45)]*(sign(D(1,1))*D(:,1)); 
        Gang = atan2d(evec1(2),evec1(1)); 
        
        % print data to output files
        fprintf(fileID1,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',timestep,sigmax2,sigmay2,sigmaxy,front_cov,pop_load,L(2,2),L(1,1),Gang);
        
        for i=1:size(pop,1)
            for j=1:size(pop,2)
                if(pop(i,j)>0)
                    fprintf(fileID2,'[%i,%i],',fitx(i),fity(j));
                    fprintf(fileID3,'%f,',pop(i,j));
                end
            end
        end
        
        fprintf(fileID2,'\n');
        fprintf(fileID3,'\n');
    end

end

v = meanfitness/steps;
v1 = meanfitx/steps;
v2 = meanfity/steps;

% close output files
fclose(fileID1);
fclose(fileID2);
fclose(fileID3);

end

