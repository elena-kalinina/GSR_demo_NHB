function [sign, u_hat]=learning_LS_new6(H, G, tr_seq, tr_length, init, y, param)

%H is the horizon - in trials (currently 5 trials each 12 secs)
%tau is the time step or time window
%x is the observed signal
%y is the desired level of signal that I currently define as the signal
%predicted with modelling the subject's negative vs neutral responses (how can we do it otherwise??)
sr=50;
%param=cfg.param;
%generate all possible sequences for a given horizon - however, we can just
%generate sequences for the number of trials e.g. 5 and then "inflate them"
%with zeros or ones according to the trial, we do not need to process all
%random sequences of 0 or 1, which reduces computations

newcontrol = generate_seq(H, [0, 1], tr_seq);  %dec2bin(0:(2^H*11)-1);
%newcontrol=[repmat(x(length(x)-tau), size(newcontrol, 1), 1) newcontrol];
%then we learn our matrix G from the first H trials of the session using
%the observed signal and stimulation protocol=control
newcontrol=newcontrol.* param;  %abs(param);
%size(newcontrol)

%newcontrol=newcontrol(:, 1:n_of_ts_b);
newcontrol=[repmat(init, size(newcontrol, 1), 1) newcontrol];

%newcontrol=[repmat([init zeros(1, 2)], size(newcontrol, 1), 1) newcontrol]; %2 number of secs before the image onset

%G=learn_G(n_of_ts_b,n_of_ts_s,tau, tau_bar, x, u);
%G_small=G(1:length(y), 1:length(y));
%G_small=G(1:size(newcontrol, 2), 1:size(newcontrol, 2));
%size(u)
newcontrol=newcontrol(:, 1:size(G, 1));
%next from the generated sequences we choose the sequence that minimizes
%the error between G*newcontrol and the response predicted by the model for
%the next 5 trials of the session (y). we know the actual protocol, so we can
%compare the learned control to the actual one 
err=zeros(size(newcontrol));
for i=1:size(newcontrol, 1)
%    figure;
%    sign=G\newcontrol(i, :)';
%    sign=resample(sign, sr, 1);
%    plot(sign)
%    hold on
%    plot([1 length(sign)],[y y])    %([1 length(newcontrol(i, :))],[y y])
% %   err(i)=sum((G_small)/newcontrol(i, :) - y); %for inv(G)
    
  err(i, :)=((G\newcontrol(i, :)')-y).^2; %for inv(G)
    
end
%err=zeros(1, size(newcontrol, 1));

% figure;
% 
% for c=1:(size(err, 1))
%     hold on
% plot(err(c, :))
% 
% end
% title('Error at different timestep of the control sequences')
% xlabel('Timestep')
% ylabel('Error')


for k=1:size(err, 1)
    new_err(k)=abs(sum(err(k, :)));
end

% figure;
% plot(new_err)
% title('Cost of different control sequences')
% xlabel('Control sequence Number')
% ylabel('Summed error of the sequence')

seqs=1:14;
%i=randint(1,1,size(newcontrol, 1)-1)+1;
for i=1:size(newcontrol, 1)
    figure;
subplot(1,2,1)
sign=G\newcontrol(i, :)';
sign=resample(sign, sr, 1);

plot(sign)
hold on
plot([1 length(sign)],[y y])    %([1 length(newcontrol(i, :))],[y y])
%   err(i)=sum((G_small)/newcontrol(i, :) - y); %for inv(G)
title('Predicted signal with respect to target signal')
xlabel('Time')
ylabel('Conductance, S(iemence)')
legend('Predicted signal', 'Target signal')
subplot(1,2,2)
plot(seqs, new_err,'o','MarkerEdgeColor','k','MarkerFaceColor','r', 'MarkerSize', 5)
set(gca,'XLim',[0 length(new_err)+1])
set(gca, 'xtick', 1:length(new_err))

%plot(new_err)
%set(gca, 'xtick', 1:length(new_err))
title('Cost of different control sequences')
xlabel('Control sequence Number')
ylabel('Summed error of the sequence')
end
% figure;
% counter=1;
% for l=1:10:length(new_err)
%     
%    if l+1<length(new_err)
%    new_err1(counter)=sum(new_err(l:l+1))
%    counter=counter+1
%    end
% end
% trials={'trial 1', 'trial 2', 'trial 3', 'trial 4', 'trial 5', 'trial 6'};
% plot(new_err1);
% title('Cost of control at different trials within the timewindow')
% set(gca,'XLim',[0 7],'XTick',1:6,'XTickLabel', trials)
%newcontrol=newcontrol(err>0, :);
%err=err(err>0);

%for trials
%newcontrol=newcontrol(2:end, :);
%err=err(2:end);
%%%%%%%%

%err=abs(sum(err, 2));
%%%%%%%%
%u_hat=find(new_err1==min(new_err1)); 

% % % err=abs(sum(err, 2));
u_hat_temp=newcontrol(new_err==min(new_err), :);

if size(u_hat_temp, 1)>1
    u_hat=u_hat_temp(randi(size(u_hat_temp, 1)), :)
    u_hat=u_hat(2:end);
else
    u_hat=u_hat_temp(2+2:length(newcontrol(new_err==min(new_err), :)));
end
%u_hat(1)=0;
u_hat=find(u_hat>0) %find(u_hat==1);

%Only for all timesteps
u_hat=unique(round(u_hat/tr_length))+1;
return



function c_seq=generate_seq(N, V, tr_seq)

% states = vec_st;
% n=length(states);


nv = length(V) ;
C = zeros(nv^N,N) ; % declaration
for ii=1:N,
     cc = 1 ;
    for jj=1:(nv^(ii-1)),
        for kk=1:nv,
            for mm=1:(nv^(N-ii)),
                C(cc,ii) = V(kk) ;
                cc = cc + 1 ;
            end
        end
    end
end


% v = mod(perms(1:(2*(m-1))),n) + 1;

%%%%%%%%%%
seq = unique(C(:,1:N),'rows');
seq(all(seq==0,2),:)=[];
seq(all(seq==1,2),:)=[];
%%%%%%%%%%%%
% seq = states(ind);

%%%%%%%%%%%%%%%%%%%
% c_seq = unique(C(:,1:N),'rows');
% c_seq(all(c_seq==0,2),:)=[];
%%%%%%%%%%%%%%%%%%%
%Only if we take the whole signal not in steps
c_seq=expand_control1(seq, tr_seq);
    

%%c_seq=expand_control(seq, tr_length);
return
    
%length of the chain
%n number of states in the chain


