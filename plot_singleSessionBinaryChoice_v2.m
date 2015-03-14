%plot_singleSessionBinaryChoice
% plots rat's data for one session


f_h = figure('Position', [556 398 745 833],'Color','w');
ax_txt = axes('position',[0.05,.87,0.4,0.12]);
set(ax_txt,'XLim',[0 1],'Ylim',[0 7]);
axis off

txt_cell = [fieldnames(session.info),struct2cell(session.info)];
%txt_cell(1:6,:); % lazy wat to get rid of 'response' field
hpos =1;
for k=1:6
    text(0.5, 7-k,[txt_cell{k,1},': '],'HorizontalAlignment','Right','FontSize', 8);
    h=text(0.5, 7-k,[txt_cell{k,2}],'HorizontalAlignment','Left','FontSize', 8);
end


R_h = axes('Position',[0.12, 0.75, 0.8, 0.10],'XTick',[],'XTickLabel',[],'Color','none','XColor','w', 'NextPlot', 'add');
m_h = axes('Position',[0.12, 0.7, 0.8, 0.05],'Visible','off','Color','none', 'NextPlot', 'add');
L_h = axes('Position',[0.12, 0.6, 0.8, 0.10],'YDir','reverse','Color','none', 'NextPlot', 'add');

%m_h = axes('Position',[0.15, 0.45, 0.8, 0.4]); % set main axes

tr = session.tr;
%plt_offset = 0; % offset for plotting - so small delay diff are discrnable

for k=1:length(tr)-1
    plot(R_h, [k-0.5,k+0.5], [tr(k).params.R_reward_delay, tr(k).params.R_reward_delay], 'y', 'LineWidth', tr(k).params.R_reward_num*4); hold on
    plot(L_h, [k-0.5,k+0.5], [tr(k).params.L_reward_delay, tr(k).params.L_reward_delay], 'y', 'LineWidth', tr(k).params.L_reward_num*4); hold on
    
    if(isempty(tr(k).well_id))
        plot(k,0,'x','r')
    else
        if(strcmpi(tr(k).well_id, 'R'))
            if(tr(k).error_flg)
                if( tr(k).well_early_go_flg)
                    plot(R_h, k, tr(k).params.R_reward_delay, 'gx', 'MarkerSize', 12);
                else
                    plot(R_h, k, tr(k).params.R_reward_delay, 'rx', 'MarkerSize', 12); % plot 'x' at expected delay if rat went right and did not get a reward 
                end
            else
                plot(R_h, k, (tr(k).R_reward_ts(1)-tr(k).well_on_ts),'ko', 'MarkerSize', 1+length(tr(k).R_reward_ts)*6); % plot 'o' at delay of reward, and scale size to number of rewards
            end
        end
        
        if(strcmpi(tr(k).well_id, 'L'))
            if(tr(k).error_flg)
                if( tr(k).well_early_go_flg)
                    plot(L_h, k, tr(k).params.L_reward_delay, 'gx', 'MarkerSize', 12);
                else
                    plot(L_h, k, tr(k).params.L_reward_delay, 'rx', 'MarkerSize', 12); % plot 'x' at expected delay if rat went right and did not get a reward 
                end
            else
                plot(L_h, k, (tr(k).L_reward_ts(1)-tr(k).well_on_ts), 'ko', 'MarkerSize', 1+length(tr(k).L_reward_ts)*6); % plot 'o' at delay of reward, and scale size to number of rewards
            end
        end
        
    end
    %odor_port_early_go_flg
    
end
linkaxes([R_h,m_h,L_h],'xy')
xlabel(L_h, sprintf('trial number'))
ylabel(L_h, sprintf('L\nReward delay (ms)'))
ylabel(R_h, sprintf('R\nReward delay (ms)'))
%ylabel(m_h, 'N.R.')
%% plot behaviour
% text: # total, # rewards, %reward, # early go, # well wait err, #R, #L,, avg lick L, avg lick R, 


%tr(k).well_on_ts-tr(k).odor_port_off_ts
    


