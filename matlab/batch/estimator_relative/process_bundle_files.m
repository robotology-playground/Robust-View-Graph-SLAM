
%close;
clc;

folder = '/home/tabuhashim/Documents/data/blue/run_20160317';
load(strcat(folder,'/options'));
filename = [folder '/bundle_5.txt'];
ncams = 20;
%color = {'r','g','b','k','c','m','y','r','g','b'};
color = hsv(ncams);
%clf;
%while(1);
    
    try a = load(filename);
        
        %clf;
        
        % addition step marker
        data = a(:, ncams*6+2);
        add = find(diff(data)>0)';
        
        % tx
        figure;%(1); clf;
        %subplot(1,3,1);
        hold on; title('tx');
        A = Inf; B = -Inf;
        for i = 1:2:ncams
            data = a(:, 6*(i-1)+1)*1000;%tx
            A = min(A, min(data)); B = max(B, max(data));
            %plot(data,color{(i+1)/2},'linewidth',2);
            plot(data,'color',color(i,:),'linewidth',2);
        end
        for i = 2:2:ncams
            data = a(:, 6*(i-1)+1)*1000;%tx
            A = min(A, min(data)); B = max(B, max(data));
            %plot(data,[color{i/2},'--'],'linewidth',2);
            plot(data,'color',color(i,:),'linewidth',2);
        end
        plot([add;add],[A*ones(1,length(add));B*ones(1,length(add))],'r');
        grid on; box on;
        
        % % ty
        % figure(2); hold on; title('ty');
        % A = Inf; B = -Inf;
        % for i = 1:2:ncams
        %     data = a(:, 6*(i-1)+2)*1000;%ty
        %     A = min(A, min(data)); B = max(B, max(data));
        %     plot(data,color{(i+1)/2},'linewidth',2);
        % end
        % for i = 2:2:ncams
        %     data = a(:, 6*(i-1)+2)*1000;%ty
        %     A = min(A, min(data)); B = max(B, max(data));
        %     plot(data,[color{i/2},'--'],'linewidth',2);
        % end
        % plot([add;add],[A*ones(1,length(add));B*ones(1,length(add))],'r');
        % grid on; box on;
        
        % % tz
        % figure(3); hold on; title('tz');
        % A = Inf; B = -Inf;
        % for i = 1:2:ncams
        %     data = a(:, 6*(i-1)+3)*1000;%tz
        %     A = min(A, min(data)); B = max(B, max(data));
        %     plot(data,color{(i+1)/2},'linewidth',2);
        % end
        % for i = 2:2:ncams
        %     data = a(:, 6*(i-1)+3)*1000;%tz
        %     A = min(A, min(data)); B = max(B, max(data));
        %     plot(data,[color{i/2},'--'],'linewidth',2);
        % end
        % plot([add;add],[A*ones(1,length(add));B*ones(1,length(add))],'r');
        % grid on; box on;
        
        % ay
        figure;%(4); clf;
        %subplot(1,3,2);
        hold on;  title('ay');
        A = Inf; B = -Inf;
        for i = 1:2:ncams
            data = a(:, 6*(i-1)+5)*180/pi;%tx
            A = min(A, min(data)); B = max(B, max(data));
            %plot(data,color{(i+1)/2},'linewidth',2);
            plot(data,'color',color(i,:),'linewidth',2);
        end
        for i = 2:2:ncams
            data = a(:, 6*(i-1)+5)*180/pi;%tx
            A = min(A, min(data)); B = max(B, max(data));
            %plot(data,[color{i/2},'--'],'linewidth',2);
            plot(data,'color',color(i,:),'linewidth',2);
        end
        plot([add;add],[A*ones(1,length(add));B*ones(1,length(add))],'r');
        grid on; box on;
        
        % % ax
        % figure(5); hold on; title('ax');
        % A = Inf; B = -Inf;
        % for i = 1:2:ncams
        %     data = a(:, 6*(i-1)+4)*180/pi;%ty
        %     A = min(A, min(data)); B = max(B, max(data));
        %     plot(data,color{(i+1)/2},'linewidth',2);
        % end
        % for i = 2:2:ncams
        %     data = a(:, 6*(i-1)+4)*180/pi;%ty
        %     A = min(A, min(data)); B = max(B, max(data));
        %     plot(data,[color{i/2},'--'],'linewidth',2);
        % end
        % plot([add;add],[A*ones(1,length(add));B*ones(1,length(add))],'r');
        % grid on; box on;
        
        % % az
        % figure(6); hold on;  title('az');
        % A = Inf; B = -Inf;
        % for i = 1:2:ncams
        %     data = a(:, 6*(i-1)+6)*180/pi;%tz
        %     A = min(A, min(data)); B = max(B, max(data));
        %     plot(data,color{(i+1)/2},'linewidth',2);
        % end
        % for i = 2:2:ncams
        %     data = a(:, 6*(i-1)+6)*180/pi;%tz
        %     A = min(A, min(data)); B = max(B, max(data));
        %     plot(data,[color{i/2},'--'],'linewidth',2);
        % end
        % plot([add;add],[A*ones(1,length(add));B*ones(1,length(add))],'r');
        % grid on; box on;
        
        figure;%(7);  clf;
        %subplot(1,3,3);
        hold on; title('number of inlier constraints');
        A = Inf; B = -Inf;
        data = a(:, end);
        A = min(A, min(data)); B = max(B, max(data));
        plot(data,'linewidth',2);
        plot([add;add],[A*ones(1,length(add));B*ones(1,length(add))],'r');
        grid on; box on;
        
        pause(10);
        %break;
        
    catch
        
        disp('Loading failed, trying again ...');
        
    end
    
%end