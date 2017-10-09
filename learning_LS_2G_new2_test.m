function [G_small, u_hat]=learning_LS_2G_new2_test(H1,n_of_timesteps_b, n_of_timesteps_s, tr_length, tau, tau_bar, x, y, u, param)

%param_a=1;

newcontrol = generate_seq(H1, [0, 1], tr_length);  %dec2bin(0:(2^H*11)-1);
newcontrol=[repmat(x(length(x)), size(newcontrol, 1), 1) newcontrol];
%[repmat(x(length(x)-tau), size(newcontrol, 1), 1) newcontrol];
%then we learn our matrix G from the first H trials of the session using
%the observed signal and stimulation protocol=control
newcontrol=newcontrol.*param;
G=learn_G(n_of_timesteps_b, n_of_timesteps_s, tau, tau_bar, x, u);
%G_small=G(1:length(y), 1:length(y));
G_small=G(1:size(newcontrol, 2), 1:size(newcontrol, 2));


%err=zeros(1, size(newcontrol, 1));
err=zeros(size(newcontrol));
for i=1:size(newcontrol, 1)
    
    %%%PLOT
%    figure;
%    plot(G_small\newcontrol(i, :)')
%    hold on
%    plot([1 length(newcontrol(i, :))],[y y])
%   err(i)=sum((G_small)/newcontrol(i, :) - y); %for inv(G)
    
  err(i, :)=((G_small\newcontrol(i, :)')-y).^2; %for inv(G)
    
end


% % % for k=1:size(err, 2)
% % %     new_err(k)=abs(sum(err(:, k)));
% % % end
% % % figure;
% % % counter=1;
% % % for l=1:24:length(new_err)
% % %     l
% % %    if l+1<length(new_err)
% % %    new_err1(counter)=sum(new_err(l:l+1))
% % %    counter=counter+1
% % %    end
% % % end
% % % trials={'trial 1', 'trial 2', 'trial 3', 'trial 4', 'trial 5'};
% % % plot(new_err1);
% % % title('Cost of control at different trials within the timewindow')
% % % set(gca,'XLim',[0 6],'XTick',1:5,'XTickLabel', trials)
% % % %newcontrol=newcontrol(err>0, :);
% % % %err=err(err>0);
% % % 
% % % %for trials
% % % %newcontrol=newcontrol(2:end, :);
% % % %err=err(2:end);
% % % 
% % % %err=err(:, 1:trial_length:size(err, 2));
% % % 
% % % err=abs(sum(err, 2));
% % % %%%%%%%%
% % % u_hat=find(new_err1==min(new_err1));    %newcontrol(err==min(err), :);
% % % %u_hat=u_hat(1:trial_length:length(u_hat));
% % % %u_hat=find(u_hat==1);
% % % %%%%u_hat=u_hat(2:end); FINDME
% % % %Only for all timesteps
% % % %%%%%u_hat=unique(round(u_hat/trial_length)+1); FINDME

%%%PLOT ERROR
% figure;
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
% set(gca, 'xtick', 1:length(new_err))
% title('Cost of different control sequences')
% xlabel('Control sequence Number')
% ylabel('Summed error of the sequence')

% figure;
% i=randint(1,1,size(newcontrol, 1)-1)+1;
% subplot(1,2,1)
% plot(G_small\newcontrol(i, :)')
% hold on
% plot([1 length(newcontrol(i, :))],[y y])
% title('Predicted signal with respect to target signal')
% xlabel('Time')
% ylabel('Conductance, S(iemence)')
% legend('Predicted signal', 'Target signal')
% subplot(1,2,2)
% plot(new_err)
% set(gca, 'xtick', 1:length(new_err))
% title('Cost of different control sequences')
% xlabel('Control sequence Number')
% ylabel('Summed error of the sequence')


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

err=abs(sum(err, 2));
%%%%%%%%
%u_hat=find(new_err1==min(new_err1)); 

% % % err=abs(sum(err, 2));
u_hat=newcontrol(err==min(err), :);
u_hat=u_hat(2:length(newcontrol(err==min(err), :)));   %u_hat(2:end);
u_hat=find(u_hat>0);  %(u_hat==1);

%Only for all timesteps
u_hat=unique(round(u_hat/tr_length))+1;
return


function G=learn_G(n_of_ts_b, n_of_ts_s, tau, tau_bar, x, u)


r=zeros(n_of_ts_b, 1);  %(round(timings(length(timings))/tau), 1); %number of timewindows ??

%define r
counterr=1;
for j_tau=1:tau:n_of_ts_s %check. source of error - mix up t counter and tau counter
    %+1
    % t
    
    for t=1:length(r)
        if t*tau_bar - j_tau >0
            r_tau(t)=u(t)*x(t*tau-j_tau);
          %  r_tau(t)=(u(t-1)*x(t*tau-j_tau))/(n_of_ts_b-t);
            %r_tau(t)=u(timings(t)-tau*j_tau+1)*x(timings(t));
        else
            %  r_tau(t)=0;
            continue;
        end
    end
    r(counterr)=sum(r_tau);
 
    counterr=counterr+1;
end

%define q
q= zeros(n_of_ts_b, 1); %zeros(round(timings(length(timings))/tau), 1); %number of timewindows ??
counter=1;
for k_tau=1:tau:n_of_ts_s
    
    for t=1:length(q)
        if t*tau_bar - k_tau >0
            q_tau(t)=x(t*tau)*x(t*tau-k_tau);
          %  q_tau(t)=(x(t*tau)*x(t*tau-k_tau))/n_of_ts_b-t;
        else
            
            %   q_tau(t)=0;
            continue;
        end
    end
    q(counter)=sum(q_tau);
   
  
    counter=counter+1;
end

%define Q
Q=zeros(n_of_ts_b, n_of_ts_s- tau_bar);  %(tau-tau_bar));
%(round(timings(length(timings))/tau));
counterq=1;
%round(timings(length(timings))/tau)
% n_tau%length(T/tau)
for n_tau_bar=1:tau:n_of_ts_s
    for n_tau=1:n_of_ts_s  %&round(timings(length(timings))/tau) %length(T/tau)
        % n_tau_bar
        for t=1:n_of_ts_b
            %  t
            if t*tau-n_tau_bar >0 && t*tau - n_tau >0
                Q_t(t, n_tau)=x(t*tau-n_tau_bar)*x(t*tau - n_tau);
            else
                % Q_t(t, n_tau_bar)=0;
                continue;
            end
            
        end
        
    end
    
    Q(counterq, :)=sum(Q_t, 1);
    counterq=counterq+1;
end

%here we get our linear system of order L - a vector
a1=Q'/(q'-r'); %in fact, it is inv(Q)*q-r
a2=a1(1:tau:length(a1));
%a2=a2*param_a;
% figure;
% plot(a2)
% title('a observed')

%%%PLOT
% figure;
% plot(a2)
% title('a observed')
%a1=(a1-min(a1))/(max(a1)-min(a1));
%here we transform it in the G matrix which is kind of
%"pseudoautoregressive" :D
%FINDME !!!!!!!!
G=autoregress(length(a2), a2); % 10 + *248); %31 7 + 3); %60 8 + 3); %245 10 + 3); %250); %159);

return


function A=autoregress(T, vec)
vec=vec*-1;
A=eye(T);

% put a(1) on -1 diagonal, a(2) on -2 etc
%n_diag=size(A, 2)-1;

for l=2:size(A, 2)
    
    if l<=length(vec)
    newdiag=repmat(vec(l-1), length(diag(A, -(l-1))), 1);
    end
    
    A=A+diag(newdiag, -(l-1));
    
end

return


function c_seq=generate_seq(N, V, tr_l)

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
%seq(all(seq==1,2),:)=[];
%%%%%%%%%%%%
% seq = states(ind);

%%%%%%%%%%%%%%%%%%%
% c_seq = unique(C(:,1:N),'rows');
% c_seq(all(c_seq==0,2),:)=[];
%%%%%%%%%%%%%%%%%%%
%Only if we take the whole signal not in steps
%%%%c_seq=expand_control(seq, tr_l); % If we are not using real onsets
    

c_seq=expand_control(seq, tr_l);

return
    
%length of the chain
%n number of states in the chain

