package beast.app.beauti;

import java.awt.Dimension;
import java.util.List;

import javax.swing.JTextField;

import beast.app.draw.InputEditor;
import beast.app.draw.ListInputEditor;
import beast.app.draw.StringInputEditor;
import beast.core.Input;
import beast.core.Logger;
import beast.core.BEASTObject;
import beast.core.BEASTInterface;



public class LoggerListInputEditor extends ListInputEditor {
    private static final long serialVersionUID = 1L;

	public LoggerListInputEditor(BeautiDoc doc) {
		super(doc);
	}

    @Override
    public Class<?> type() {
        return List.class;
    }

    @Override
    public Class<?> baseType() {
        return Logger.class;
    }
    

    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
    	super.init(input, plugin, itemNr, bExpandOption, bAddButtons);
    }
    
    @Override
    protected void addSingleItem(BEASTInterface plugin) {
    	currentLogger = (Logger) plugin;
    	super.addSingleItem(plugin);
    }
    
    public InputEditor createFileNameEditor() throws Exception {
        Input<?> input = currentLogger.fileNameInput;
        StringInputEditor fileNameEditor = new StringInputEditor(doc);
        fileNameEditor.init(input, currentLogger, -1, ExpandOption.FALSE, true);

        // ensure file name entry has larger size than the standard size
        JTextField fileNameEntry = fileNameEditor.getEntry();
        Dimension size = new Dimension(400, fileNameEntry.getPreferredSize().height);
        fileNameEntry.setMinimumSize(size);
        fileNameEntry.setPreferredSize(size);
        return fileNameEditor;
    }
    
    Logger currentLogger;
}
