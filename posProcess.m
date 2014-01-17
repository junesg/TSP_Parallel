%



filename = {'lin105', 'rat195','eil101', 'a280', 'kroB200',...
    'kroC100','kroD100', 'kroE100', 'pcb442','pla7397','rd100','pr1002'};
tour_len = [105, 195, 101,280,200,100,100,100,443,7397,100,1002];
TSP_result = [14379, 2323, 629, 2579, 29437, 20749, 21294,22068, 50778, 23260728, 7910,259045];
counter = 1;



for iter = 10:1:10% length(tour_len)
    
    
    directory = strcat('output/', filename{iter},'.o*');
    output_files = dir(directory);
    
    
    for i = 1:1:length(output_files)
        output_files(i).name
        
        text = fileread(strcat('output/',output_files(i).name));
        %text = fileread('output/lin105.o11484915');%4735');
        aText = regexprep(text,'[^\w.=]','');

        ppn_exp = '(np=|ppn=|procs=|prcs=)+[\d*]+([I])';
        node_exp = 'nodes=+[\d*]+[pn]';
        distance_exp = 'distance=[\d*\.]+[\d*]';
        time_exp = 'TimeSpent[\d*\.]+[\d*]';


        fileread_info_node = regexp(aText, node_exp, 'match');
        fileread_info_ppn = regexp(aText, ppn_exp, 'match');
        fileread_info_dist = regexp(aText, distance_exp, 'match');
        fileread_info_time = regexp(aText, time_exp, 'match');


        if (isempty(fileread_info_ppn))
            output_files(i).name;
            error('processor missing');
        end
        %if the file is a non success, do nothing
        %if the file is not empty, then we store
        if (~isempty(fileread_info_dist) || ~isempty(fileread_info_time))
            
            distance = zeros(1,length(fileread_info_dist));
            
            for index = 1:1:length(fileread_info_dist)
                distText = regexprep(fileread_info_dist{index},'[a-z=]','');
                distance(index) = str2double(distText);
            end
            
            time = str2double(regexprep(fileread_info_time,'[a-z=A-Z]',''));
            if(isempty(fileread_info_node))
                node = 1;
            else
                node = str2double(regexprep(fileread_info_node,'[a-z=A-Z]',''));
            end
            procs = str2double(regexprep(fileread_info_ppn,'[a-z=A-Z]',''));
            processors = node * procs;
            accuracy = 1- (min(distance)-TSP_result(iter))/TSP_result(iter);
            
            data{counter}{1} = filename{iter}; 
            data{counter}{2} = tour_len(iter);
            data{counter}{3} = processors;
            data{counter}{4} = time;
            data{counter}{5} = max(distance);
            data{counter}{6} = min(distance);
            data{counter}{7} = mean(distance);
            data{counter}{8} = accuracy;
            data{counter}{9} = TSP_result(iter);
            counter = counter +1;
            
            fid = fopen('Output.txt','a');
            if fid~=-1
                fprintf(fid, '%d %s %d %f %f %f %f %f %f %f\n',iter,filename{iter},tour_len(iter),...
                    processors, time, max(distance), min(distance), mean(distance), accuracy, TSP_result(iter));
            else
                error('unopened');
            end
            fclose(fid);
        end


    end
end













