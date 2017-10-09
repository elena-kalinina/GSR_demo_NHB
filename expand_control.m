function new_seq=expand_control(seq, tr_l)

for s=1:size(seq, 1)
    
    c_seq_temp=[];
    for tr=1:length(seq(s, :))
        %    tr
        myvec=zeros(1, tr_l);
        if seq(s, tr)==1
            %   myvec=ones(1, tr_l);
            %myvec(1:round(tr_l/3))=1;
            myvec(1)=1;
            % else
            %    myvec=zeros(1, tr_l);
        end
        % c_seq_temp=[c_seq_temp seq(s, tr) myvec];
        %  length(myvec)
        c_seq_temp=[c_seq_temp myvec];
        
    end
    new_seq(s, :)=c_seq_temp;
    
end



return


