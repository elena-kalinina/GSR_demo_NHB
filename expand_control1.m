function new_seq=expand_control1(seq, tr_seq)

for s=1:size(seq, 1)
    tr_seq;
    c_seq_temp=zeros(1, tr_seq(end));
    length(c_seq_temp);
   % for tr=1:size(seq, 1)
        %    tr
     %   myvec=zeros(1, tr_seq);
     
     %   if seq(s, tr)==1
            
            ind=find(seq(s, :)==1);
            %   myvec=ones(1, tr_l);
            %myvec(1:round(tr_l/3))=1;
            c_seq_temp(tr_seq(ind))=1;
            % else
            %    myvec=zeros(1, tr_l);
      %  end
        % c_seq_temp=[c_seq_temp seq(s, tr) myvec];
        %  length(myvec)
      
        
  %  end
    new_seq(s, :)=c_seq_temp;
    
end



return


