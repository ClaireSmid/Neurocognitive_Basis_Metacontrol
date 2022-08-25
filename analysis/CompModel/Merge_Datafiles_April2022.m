
load kidsresults_P6_T0_5April2022.mat

Final_out = [];

Final_out(:).ID = int32([kidresults.id]');
%Final_out(:).nstarts = [kidresults.nstarts]';

% extract P6 w overall
a = cell2mat({kidresults.x}.');
Final_out.it_P6 = a(:,1);
Final_out.lr_P6 = a(:,2);
Final_out.eg_P6 = a(:,3);
Final_out.w_P6 = a(:,4);
Final_out.st_P6 = a(:,5);
Final_out.repst_P6 = a(:,6);

load kidsresults_P7_T0_5April2022.mat

% extract P6 w overall
b = cell2mat({kidresults.x}.');
Final_out.it = b(:,1);
Final_out.lr = b(:,2);
Final_out.eg = b(:,3);
Final_out.w_lo = b(:,4);
Final_out.w_hi = b(:,5);
Final_out.st = b(:,6);
Final_out.repst = b(:,7);

% save output file
T = struct2table(Final_out);

writetable(T,'MB_MF_P6_P7_11April2022.csv','Delimiter',',')