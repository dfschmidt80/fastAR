% Black and white or Colour?
if (~exist('BW', 'var'))
    BW = false;
end

% Large (single) figure properties
LARGE_FIG_HEIGHT = 3.5;
LARGE_FIG_WIDTH  = 6.5;
LARGE_FIG_AXIS_FONTSIZE = 11;
LARGE_FIG_LEGEND_FONTSIZE = 10;
LARGE_FIG_FONTNAME = 'Times New Roman';
LARGE_FIG_LINEWIDTH = 1.5;
LARGE_FIG_MARKERSIZE = 6;

LARGE_THREE_FIG_HEIGHT = 2.5;

% Small (multi) figure properties
SMALL_FIG_HEIGHT = 5.8*0.45;
SMALL_FIG_WIDTH  = 6.8*0.45;
SMALL_FIG_LINEWIDTH = 0.9;
SMALL_FIG_AXIS_FONTSIZE = 10;
SMALL_FIG_LEGEND_FONTSIZE = 9;
SMALL_FIG_FONTNAME = 'Times New Roman';
SMALL_FIG_MARKERSIZE = 5;

% Colour definitions
% Black and white
if (BW)
    %myblack = [0 0 0];
    %mydarkgrey = [0.4, 0.4, 0.4];
    %mygrey  = [0.65 0.65 0.65];
    
    myblue = [0 0 0];
    mybluestyle = '-';
    mybluediscretestyle = '-o';
    
    myorange = [0 0 0];
    myorangestyle = '--';
    myorangediscretestyle = '--o';
    
    myyellow = [0.65, 0.65, 0.65];
    myyellowstyle = '-';
    myyellowdiscretestyle = '-o';
    
    mygreen = [0.65, 0.65, 0.65];
    mygreenstyle = '--';
    mygreendiscretestyle = '--o';

% Colour
else
    myblue = [0 114 189]/256;
    mybluestyle = '-';
    mybluediscretestyle = '-o';

    myorange = [0.85, 0.33, 0.1];
    myorangestyle = '-';
    myorangediscretestyle = '-o';
    
    myyellow = [237 177 32]/256;
    myyellowstyle = '-';
    myyellowdiscretestyle = '-o';
    
    mygreen = [119 172 48]/256;
    mygreenstyle = '-';
    mygreendiscretestyle = '-o';

    mypurple = [126 47 142]/256;
    mypurplestyle = '-';
    mypurplediscretestyle = '-o';
end