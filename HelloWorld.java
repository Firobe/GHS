import jbotsim.Topology;
import jbotsim.ui.JViewer;

public class HelloWorld{
    public static void main(String[] args){
        Topology tp = new Topology(false);
        tp.setDefaultNodeModel(BasicNode.class);
        Import.importGraph(tp, "pouet.dot", 150);
        tp.start();
    }
}    
