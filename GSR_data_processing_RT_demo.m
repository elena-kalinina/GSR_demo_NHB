function void = GSR_data_processing_RT_demo(SubjID, ExpType, sesN, cfg)
H1=4; %6;
H2=3; %4;
n_trial=H2;
predicted_signal=[];
sr=10; %cfg.sr;
cfg.filter.sr=49; %50; %51.2;
counter=1;
tau=cfg.tau;
n_tau=n_trial*cfg.trial_length;
n_tau_bar=n_trial*cfg.trial_length*tau;
tau_bar=cfg.tau_bar;
bl_dur=cfg.bl_dur;
async=0; %cfg.async;
tr_length=cfg.trial_length;
pr_filename=sprintf('%s_%s-%d.mat', SubjID, ExpType, sesN);
protocol=load(pr_filename, 'onsets', 'control');
protocol.onsets;
filename=sprintf('%s%d.mat', SubjID, sesN);
load(filename, 'data');
mydata1=data{1};
new_onsets1= (protocol.onsets-bl_dur);
new_onsets2=protocol.onsets*sr+1;
new_onsets3=round((protocol.onsets-bl_dur)*cfg.filter.sr)+1; %
new_onsets_red=[new_onsets1(1); new_onsets1(H1+1:H2:length(new_onsets1))]+4;
f=1;
data_start=round(cfg.to_cut*cfg.filter.sr)+1;

predicted_control=zeros(length(new_onsets1), 1);

%%

while 1
    
    mystart=input('enter when you are ready to start \n');
    if mystart==1
        
        
        
        timepoint=120;
        
        alert=0;
        elapsedTime = 0;
        tic;% Reset to 0
        session_start=tic;
        mytime=0;% TRUE if the shimmer starts streaming
        
        plotData = [];
        timeStamp = [];
        mydata_temp=[];
        h.figure1=figure('Name','Shimmer 1 signals');
        set(h.figure1, 'Position', [100, 500, 800, 400]);
        
        while timepoint <=350
            
            if timepoint<1
                fprintf('No new data \n')
            else% TRUE if new data has arrived
                pause(1);
                
                length(mydata1);
                
                
                set(0,'CurrentFigure',h.figure1);
                
                plot(1:timepoint*cfg.filter.sr,mydata1(1:timepoint*cfg.filter.sr),'b')
             %   myfilm(counter)=getframe;
                counter=counter+1;
                
                if  ismember (timepoint, new_onsets_red)
                    
                    for i=1:length(predicted_control)
                        if  (protocol.onsets(i)*cfg.filter.sr)<=(timepoint+40)*cfg.filter.sr
                            if predicted_control(i)==1
                                
                                p1=line([protocol.onsets(i)*cfg.filter.sr protocol.onsets(i)*cfg.filter.sr],[0.0001 0.005],  'Color', 'm', 'LineStyle', ':', 'LineWidth', 2.2);
                                
                            else
                                p2=line([protocol.onsets(i)*cfg.filter.sr protocol.onsets(i)*cfg.filter.sr],[0.0001 0.005],  'Color', 'c', 'LineStyle', ':', 'LineWidth', 2.2);
                                
                                hold on
                            end
                        end
                    end
                  %  hold off
                    legend([p1, p2], 'negative stim', 'neutral stim', 'AutoUpdate','off')
                    xlabel('Time; sampling rate  = 50 Hz')
                    ylabel('Conductance (Siemens)')
                    
                  %  myfilm(counter)=getframe;
                    counter=counter+1;
                end
                
                
                timepoint=timepoint+1
                
                if  ismember (timepoint, new_onsets_red)
                    
                    length(mydata1);
                    mydata2=mydata1(data_start*1.5:timepoint*51);
                    figure;
                    
                    
                    plot(mydata2)
                    title('Data processed Moving Average')
                   % fprintf('eccoci\n')
                    %hold off
                    [sts, mydata, newsr] = myscr_prepdata(mydata1(data_start:timepoint*51), cfg.filter, cfg);
                    
                    figure;
                   % set(gca,'fontsize',20)
                    hold on
                    plot(mydata)
                    hold off
                    title('Data processed with Butterworth filter')
                    length(mydata);
                    
                    
                    x0=mydata(end);
                    
                    filename_o=cfg.OnsetFile;
                    load(filename_o, 'real_onsets');
                    new_onsets4=round(real_onsets(2:end))+bl_dur+async;
                    new_onsets4=new_onsets4-cfg.to_cut;
                    
                    
                    if timepoint==new_onsets_red(1)
                        baseline_mean=mean (mydata(cfg.to_cut*sr:end));
                        
                        baseline_std = std(mydata(cfg.to_cut*sr:end));
                        
                        if cfg.baseline_mean
                            predicted_resp=baseline_mean+((cfg.level+cfg.factor)*baseline_std);
                        else
                            predicted_resp=cfg.mymean+((cfg.level+cfg.factor)*cfg.mystd);
                        end
                        
                        window=1:H1+1;
                        control_seq=(new_onsets1(window)-new_onsets1(window(1)))+3; % 1 trial start 5 protocol update 3 OR 4 ??
                        load(cfg.MatrixFile);
                        G=G_small;
                        [pred_sign, u_hat]=learning_LS_new6(H1, G, control_seq, tr_length, x0, predicted_resp, cfg.param);
                        predicted_control(u_hat)=1;
                        predicted_signal=[predicted_signal pred_sign];
                    else
                        timepoint;
                        timepoint-4;
                        ind=find(new_onsets1==timepoint-4);
                        window=ind-H2:ind;
                        
                        
                        timings=new_onsets4(window)*sr+1;
                        
                        if length(mydata)>=timings(end)
                            x=mydata(timings(1)-1:timings(end)-1);
                        else
                            fprintf('Warning: data length is behind acquisition timings\n');
                            x=mydata(timings(1):length(mydata));
                        end
                        
                        
                        new_onsets_red;
                        control_seq=(new_onsets4(window)-new_onsets4(window(1)))+3;
                        timeline=zeros(length(x), 1);
                        
                        filename_c=cfg.FeedbackFile;
                        load(filename_c, 'predicted_control');
                        new_controls=predicted_control(window);
                        onsets_neg=timings(new_controls==1);
                        if ~isempty(onsets_neg)
                            onsets_neg(1);
                            onsets_neg=onsets_neg-onsets_neg(1)+1;
                            
                            for tm=1:length(timeline)
                                if ismember(tm, onsets_neg)
                                    timeline(tm)=1;
                                end
                            end
                            
                            timeline=[x(1); timeline(1:length(x)-1)];   %(1)=x(1);
                            length(x)
                            onsets_neg
                            response_amp=calculate_amplitude(x, onsets_neg, tr_length*sr);
                            [pred_sign, u_hat]=learning_LS_new4(H2, control_seq, x0, tr_length, tau, tau_bar, x, predicted_resp, timeline, response_amp);
                            
                        else
                            fprintf('amplitude error \n');
                            timeline=[x(1); timeline(1:length(x)-1)];
                            [pred_sign, u_hat]=learning_LS_new4(H2, control_seq, x0, tr_length, tau, tau_bar, x, predicted_resp, timeline, cfg.param);
                        end
                        
                        predicted_signal=[predicted_signal; pred_sign];
                        
                        test1=u_hat+H1+H2*(f-1);
                        test2=u_hat+H1+H2*(f-2);
                        H1+(H2*(f-1))+1:H1+H2*f;
                        
                        predicted_control(u_hat+H1+H2*(f-2))=1;
                        
                        
                    end
                    
                    
                    filename=cfg.FeedbackFile;
                    save(filename, 'predicted_control');
                    f=f+1;
                    clear mydata
                    
                end
                
                
                
            end
            
            elapsedTime = elapsedTime + toc;
            tic;
        end
        elapsedTime = elapsedTime + toc
        fprintf('DONE \n');
        
        break
    end
end
predicted_file=sprintf('%s_%s_pred_sign_%d.mat', SubjID, ExpType, sesN);
save(predicted_file, 'predicted_signal')

%movie(myfilm);
%end






