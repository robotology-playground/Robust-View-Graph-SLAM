close all

figure ( 9 ) ;

subplot(1,2,1);
plot(neck(1,:)'*180/pi,'k','linewidth', 2);grid on;hold on;
plot(pose1(6,1:2:end)'*180/pi,'b','linewidth', 2);
plot(pose2(6,1:2:end)'*180/pi,'g','linewidth', 2);
plot(pose3(6,1:2:end)'*180/pi,'r--','linewidth', 2);
plot(pose4(6,1:2:end)'*180/pi,'r','linewidth', 2);
legend('Tilt encoder','Kinematics','Averaging','MST','Robust');
% Add lines
%h1 = line([2 2],[1 10]);
%h2 = line([5 5],[1 10]);
% Set properties of lines
%set([h1 h2],'Color','k','LineWidth',1)
% Add a patch
patch([2 7 7 2],[-20 -20 40 40],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([15 23 23 15],[-20 -20 40 40],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([28 36 36 28],[-20 -20 40 40],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([40 46 46 40],[-20 -20 40 40],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([67 76 76 67],[-20 -20 40 40],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([78 84 84 78],[-20 -20 40 40],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([91 98 98 91],[-20 -20 40 40],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([113 118 118 113],[-20 -20 40 40],'r','FaceAlpha',0.2,'EdgeAlpha',0)
% The order of the "children" of the plot determines which one appears on top.
% I need to flip it here.
%set(gca,'children',flipud(get(gca,'children')))
axis tight

subplot(1,2,2); 
plot(stereo_pose_k(4,:)*180/pi,'b','linewidth', 2); hold on;
plot(stereo_pose_a(4,:)*180/pi,'g','linewidth', 2);
plot(stereo_pose_m(4,:)*180/pi,'m','linewidth', 2);
plot(stereo_pose_n(4,:)*180/pi,'r','linewidth', 2);
legend('Kinematics','Averaging','MST','Robust');
plot(stereo_pose_k(6,:)*180/pi,'b','linewidth', 2);
plot(stereo_pose_a(6,:)*180/pi,'g','linewidth', 2);
plot(stereo_pose_m(6,:)*180/pi,'m','linewidth', 2);
plot(stereo_pose_n(6,:)*180/pi,'r','linewidth', 2);
patch([2 7 7 2],[-1.2 -1.2 2.7 2.7],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([15 23 23 15],[-1.2 -1.2 2.7 2.7],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([28 36 36 28],[-1.2 -1.2 2.7 2.7],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([40 46 46 40],[-1.2 -1.2 2.7 2.7],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([67 76 76 67],[-1.2 -1.2 2.7 2.7],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([78 84 84 78],[-1.2 -1.2 2.7 2.7],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([91 98 98 91],[-1.2 -1.2 2.7 2.7],'r','FaceAlpha',0.2,'EdgeAlpha',0)
patch([113 118 118 113],[-1.2 -1.2 2.7 2.7],'r','FaceAlpha',0.2,'EdgeAlpha',0)
axis tight


% subplot(2,2,2); 
% plot(neck(3,:)'*180/pi,'k','linewidth', 2);grid on;hold on;
% plot(pose1(5,1:2:end)'*180/pi,'b','linewidth', 2);
% plot(pose2(5,1:2:end)'*180/pi,'g','linewidth', 2);
% plot(pose3(5,1:2:end)'*180/pi,'r--','linewidth', 2);
% plot(pose4(5,1:2:end)'*180/pi,'r','linewidth', 2);
% legend('Heading encoder','Kinematics','Averaging','MST','Robust');
% patch([3 15 15 3],[-30 -30 60 60],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([17 27 27 17],[-30 -30 60 60],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([28 38 38 28],[-30 -30 60 60],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([40 46 46 40],[-30 -30 60 60],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([53 72 72 53],[-30 -30 60 60],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([78 85 85 78],[-30 -30 60 60],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([103 108 108 103],[-30 -30 60 60],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([112 119 119 112],[-30 -30 60 60],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([121 125 125 121],[-30 -30 60 60],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% axis tight
% 
% 
% subplot(2,2,4); 
% plot(stereo_pose_k(5,:)*180/pi,'b','linewidth', 2); hold on;
% plot(stereo_pose_a(5,:)*180/pi,'g','linewidth', 2);
% plot(stereo_pose_m(5,:)*180/pi,'m','linewidth', 2);
% plot(stereo_pose_n(5,:)*180/pi,'r','linewidth', 2);
% legend('Kinematics','Averaging','MST','Robust');
% patch([3 15 15 3],[-.2 -.2 2.2 2.2],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([17 27 27 17],[-.2 -.2 2.2 2.2],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([28 38 38 28],[-.2 -.2 2.2 2.2],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([40 46 46 40],[-.2 -.2 2.2 2.2],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([53 72 72 53],[-.2 -.2 2.2 2.2],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([78 85 85 78],[-.2 -.2 2.2 2.2],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([103 108 108 103],[-.2 -.2 2.2 2.2],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([112 119 119 112],[-.2 -.2 2.2 2.2],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% patch([121 125 125 121],[-.2 -.2 2.2 2.2],'b','FaceAlpha',0.2,'EdgeAlpha',0)
% axis tight
% 
