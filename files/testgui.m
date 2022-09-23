function testgui
%GUI figure
fig=figure('Visible','off','Position',[360,500,450,285],'MenuBar','none','NumberTitle','off');

action1=uicontrol('Style','pushbutton','String','Action1','Position',[350,220,70,25],'Callback',{@button1_callback});
action2=uicontrol('Style','pushbutton','String','Action1','Position',[350,190,70,25],'Callback',{@button2_callback});
closegui=uicontrol('Style','pushbutton','String','close','Position',[350,130,70,25],'Callback',{@button3_callback});

hpopup=uicontrol('Style','popupmenu','String',{'S1','S2','S3'},'Position',[350,160,70,25],'Callback',{@popup_callback});



%align buttons along centers
align([action1,action2,hpopup],'Center','None');

%change units of components to resize automatically
set([fig,action1,action2,hpopup],'Units','normalized')
set(fig,'Name','Testgui');
movegui(fig,'center')
set(fig,'Visible','on')

%initialize popup menue
shownumber=1;
set(hpopup,'Value',3);

%POPUP-MENU CALLBACK FUNCTION

    function popup_callback(source,data)
        str=get(source,'String');
        val=get(source,'Value');
       
        switch str{val}
            case 'S1'
                shownumber=1;
            case 'S2'
                shownumber=2;
            case 'S3'
                shownumber=3;
        end
        source
        data
    end

%PUSHBUTTIONS
    function button1_callback(source,eventdata)
        fprintf(['pushed button 1. shownumber=',num2str(shownumber),'.\n']);
    end
    function button2_callback(source,eventdata)
        fprintf(['pushed button 2. shownumber=',num2str(shownumber),'.\n']);
    end
    function button3_callback(source,eventdata)
        display Goodbye
        close(fig)
    end


end