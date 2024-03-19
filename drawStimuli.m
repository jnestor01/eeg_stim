function drawStimuli(p, window, stim, destRect, fixColor, orientation, phMask, phRect, white, cx, cy)

% stim = grating tex. stim empty for no stim
% destRect = grating location rect. empty for no stim
% orientation = grating orientation

ringColor = p.ringColor*white;
bgColor = p.backgroundColor*white;

if p.showPlaceholders==1 % angular frame
    Screen('FillRect', window, bgColor);
    
    Screen('FrameRect', window, color, rect, penWidth);
    Screen('FillRect', window, bgColor, squeezeRect(rect, [0.8 1.2])) % tall rect
    Screen('FillRect', window, bgColor, squeezeRect(rect, [1.2 0.8])) % wide rect
    if ~isempty(stim)
        Screen('DrawTexture', window, stim, [], destRect, orientation);        
    end
    drawFixation(window, cx, cy, p.fixSize, fixColor*white);
end

if p.showPlaceholders==2 % ring frame
    Screen('FillRect',window, ringColor);
    if ~isempty(stim)
        Screen('DrawTexture', window, stim, [], destRect, orientation);        
    end
    Screen('DrawTexture', window, phMask)
    drawFixation(window, cx, cy, p.fixSize, fixColor*white);    
end