function display_brain_img_sub(img, background_img, z_axis, x_axis, title_name, pv, ax)
    % Displays a foreground image (img) on top of a background image.
    % Adjust the code as necessary.
    
    if nargin < 7 || isempty(ax)
        figure; 
        ax1 = axes;
    else
        axes(ax);
        ax1 = ax;
    end
    

        % Check the maximum and minimum values of img
        % max_val = max(img(:));
        % min_val = min(img(:));
        % if abs(max_val) < abs(min_val)
        %     img = -img;  % Invert the image if the maximum absolute value is less than the minimum
        % end

        % check the sign of mean value of img
        if mean(img(:)) < 0
            img = -img;  % Invert the image if the mean value is negative
        end

    % Display the background image with adjusted alpha value
    hb1 = imagesc(ax1, x_axis, z_axis, background_img); 
    colormap(ax1, 'gray'); 
    hold on;
    set(hb1, 'AlphaData', 0.4);  % Set alpha value for background image
    
    ax2 = axes; 
    hb2 = imagesc(ax2, x_axis, z_axis, img);
    
    switch pv
        case 1
            hb2.AlphaData = 5 * img; 
            cmap = 'hot';
            ax2.Visible = 'off'; 
            c = colorbar;  
            c.FontSize = 12;
            c.Ticks = linspace(min(nonzeros(img(:))), max(img(:)), 4);
            c.TickLabels = round(linspace(min(nonzeros(img(:))), max(img(:)), 4), 2);
            set(get(c, 'Title'), 'String', 'PCC');
            caxis([min(nonzeros(img(:))) max(img(:))]);
    
        case 2
            cmap = lines(1); 
            hb2.AlphaData = .4 * img;
    end
    
    colormap(ax2, cmap);
    ax2.Visible = 'off';
    ax2.YDir = 'reverse';
    linkprop([ax1, ax2], 'Position');
    % xlabel(ax1, 'Width [mm]'); 
    % ylabel(ax1, 'Depth [mm]');
    title(ax1, title_name);
    ax1.YTick = 2:7; 
    set(ax1, 'FontSize', 14);
    
    end
    
    







% function display_brain_img_sub(img, background_img, z_axis, x_axis, title_name, pv, ax)
%     % Displays a foreground image (img) on top of a background image
%     % You may adjust the code if you find it necessary
    
%     if nargin < 7 || isempty(ax)
%         figure; 
%         ax1 = axes;
%     else
%         axes(ax);
%         ax1 = ax;
%     end
    
%     imagesc(ax1, x_axis, z_axis, background_img); 
%     colormap(ax1, 'gray'); 
%     hold on;
    
%     ax2 = axes; 
%     hb = imagesc(ax2, x_axis, z_axis, img);
    
%     switch pv
%         case 1
%             hb.AlphaData = 5 * img; 
%             cmap = 'hot';
%             ax2.Visible = 'off'; 
%             c = colorbar;  
%             c.FontSize = 12;
%             c.Ticks = linspace(min(nonzeros(img(:))), max(img(:)), 4);
%             c.TickLabels = round(linspace(min(nonzeros(img(:))), max(img(:)), 4), 2);
%             % set(get(c, 'Title'), 'String', 'PCC');
%             caxis([min(nonzeros(img(:))) max(img(:))]);
    
%         case 2
%             cmap = lines(1); 
%             hb.AlphaData = .4 * img;
%     end
    
%     colormap(ax2, cmap);
%     ax2.Visible = 'off';
%     ax2.YDir = 'reverse';
%     linkprop([ax1, ax2], 'Position');
%     % xlabel(ax1, 'Width [mm]'); 
%     % ylabel(ax1, 'Depth [mm]');
%     % title(ax1, title_name);
%     ax1.YTick = 2:7; 
%     set(ax1, 'FontSize', 10);
    
%     end
    