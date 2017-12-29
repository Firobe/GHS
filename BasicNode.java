import jbotsim.*;
import jbotsim.Message;
import java.util.*;

public class BasicNode extends Node{
    public enum Status {
        INIT,
        BEGIN,
        READY,
        ROOT;
    }
    public class Patatoide {
        public Node sender;
        public Integer frag;
        public Patatoide(Node s, Integer f) {
            sender = s; frag = f;
        }
    }

    private Integer frag;
    private Node father;
    private ArrayList<Node> sons;
    private Status status;
    private Integer phase;
    private List<Node> neighbors;

    private Integer pongCounter;

    private Integer fragCounter;

    private Integer mcoeCounter;
    private Link curMCOE;
    private Node goodboy; //Son that sent minMCOE

    private ArrayList<Patatoide> mergeList;

    private Integer echoCounter;
    @Override
    public void onStart() {
        // initialize the node variables
        frag = getID();
        father = this;
        sons = new ArrayList<Node>();
        mergeList = new ArrayList<Patatoide>();
        status = Status.INIT;
        phase = 0;
        neighbors = getNeighbors();
    }

    public void init() {
        mcoeCounter = 0;
        pongCounter = 0;
        goodboy = null;
        mergeList.clear();
        echoCounter = 0;
        fragCounter = 0;
    }

    private boolean greaterThan(Link l) {
        //STRICT ORDER
        if(l == null) return false;
        if(curMCOE == null) return true;
        int comp = l.compareTo(curMCOE);
        if(comp < 0) return true;
        if(comp > 0) return false;
        int x1 = Math.min(curMCOE.endpoint(0).getID(),
                curMCOE.endpoint(1).getID()); 
        int y1 = Math.max(curMCOE.endpoint(0).getID(),
                curMCOE.endpoint(1).getID()); 
        int x2 = Math.min(l.endpoint(0).getID(),
                l.endpoint(1).getID()); 
        int y2 = Math.max(l.endpoint(0).getID(),
                l.endpoint(1).getID()); 
        return (x2 < x1 || (x1 == x2 && y2 < y1));
    }

    private void disp() {
        for(Node s : sons) {
            System.out.println(getID() + "--"
                    + s.getID() + ";");
            send(s, new Message(null, "DISP"));
        }
    }

    private void sendMCOE() {
        if(father != this)
            send(father, new Message(curMCOE, "MCOE"));
        else { //I am root
            if(curMCOE == null) {//Step 4
                System.out.println("ALGOR THERMAINT");
                System.out.println("graph {");
                disp();
            }
            else {
                //Step 5
                ackProcess(curMCOE);
            }
        }
    }

    private void mcoeProcess(Link l, Node sender) {
        //Inside step 3
        if(sender != null)
            mcoeCounter++;
        if(sender != null && greaterThan(l)) {
            curMCOE = l;
            goodboy = sender;
        }

        //Received MCOE from every son, send to father
        if(mcoeCounter == sons.size() && fragCounter == 0)
            sendMCOE();
    }

    public class MessNPhase {
        public Object object;
        public Integer phase;
        public MessNPhase(Object o, Integer p) {
            object = o;
            phase = p;
        }
    }

    @Override
    public void send(Node d, Message m) {
        MessNPhase mp = new MessNPhase(m.getContent(), phase);
        super.send(d, new Message(mp, m.getFlag()));
        System.out.println("PHASE " + phase +
                "(" + status.toString() + ")" +
                " : " + getID() + " send to " +
                d.getID() + " : " + m.getFlag());
    }

    private void nextPulse() {
        for(Node n : neighbors)
            send(n, new Message(null, "PING"));
    }

    private void processPulse() {
        //Step 1
        init();
        nextPulse();
        //processPong();
    }

    private void processPong() {
        status = Status.BEGIN;
        curMCOE = null;
        for(Node n : neighbors)
            if(sons.contains(n))
                send(n, new Message(null, "PULSE"));
            else if (n != father){
                send(n, new Message(null, "FRAG"));
                fragCounter++;
            }
    }

    private void ackProcess(Link rMCOE) {
        //Step 5
        for(Node n : sons)
            send(n, new Message(rMCOE, "ACK"));
        //Update frag number
        frag = rMCOE.endpoint(0).getID();
        status = Status.READY;
        if(curMCOE != null && rMCOE.equals(curMCOE)) { 
            //On the path to MCOE
            if(father != this) //Not current root
                sons.add(father);
            father = goodboy;
            if(goodboy != null) //If not new root, swap
                sons.remove(sons.indexOf(goodboy));
            else { //I am u*
                father = this;
                //Step 6
                send(curMCOE.endpoint(1), new Message(
                            frag, "MERGE"));
            }
        } 

        //Step 7
        for(Patatoide p : mergeList) {
            sons.add(p.sender);
            if(father == this
                && rMCOE.equals(new Link(this, p.sender))
                && frag < p.frag) {
                status = Status.ROOT;
            }
        }
        if(status == Status.ROOT) {
            System.out.println(getID() + " is now fucking root");
            //Step 8
            for(Node n : sons)
                send(n, new Message(frag, "NEW"));
        }
        curMCOE = rMCOE;
    }

    @Override
    public void onClock() {
        // code to be executed by this node in each round
        List<Message> received = getMailbox();
        //Initial pulse for the first phase
        if(status == Status.INIT) {
            processPulse();
            status = Status.BEGIN;
        }

        //Reading messages
        for(Iterator<Message> it = received.iterator() ; it.hasNext() ; ) {
            Message m = it.next();
            MessNPhase mp = (MessNPhase) m.getContent();
            if(mp.phase <= phase) {
                String flag = m.getFlag();
                Object content = mp.object;
                Node sender = m.getSender();

                //PULSE
                if(flag == "PULSE")
                    processPulse();

                else if(flag == "PING")
                    send(sender, new Message(frag, "PONG"));

                else if(flag == "PONG") {
                    pongCounter++;
                    if(pongCounter == neighbors.size())
                        processPong();
                }

                //FRAG
                else if(flag == "FRAG") //Step 2
                    send(sender,
                            new Message(frag, "RFRAG"));
                
                //RFRAG
                else if(flag == "RFRAG") { //Step 2
                    fragCounter--;
                    Link l = new Link(this, sender);
                    //Update min if edge is outgoing
                    if((Integer) content != frag
                            && greaterThan(l))
                        curMCOE = l;

                    //If already received every son MCOE and this
                    //is last RFRAG, propagate MCOE
                    if(mcoeCounter == sons.size() && fragCounter == 0)
                        sendMCOE();
                }

                //MCOE
                else if(flag == "MCOE") {
                    mcoeProcess((Link) content, sender);
                }

                //ACK
                else if(flag == "ACK") {
                    ackProcess((Link) content);
                }

                //MERGE
                else if(flag == "MERGE") {
                    if(status == Status.BEGIN ) {
                        //Step 7
                        mergeList.add(new Patatoide(sender,
                                    (Integer) content));
                    }
                    else if(status == Status.ROOT) {
                        //Step 9
                        sons.add(sender);
                        send(sender, new Message(frag, "NEW"));
                    }
                    else if(status == Status.READY) { //READY
                        Integer senderFrag = (Integer) content;
                        Link outEdge = new Link(this, sender);
                        sons.add(sender);
                        if(father == this
                            && curMCOE.equals(outEdge)
                            && frag < senderFrag) {
                            status = Status.ROOT;
                        }
                        if(status == Status.ROOT) {
                            System.out.println(getID() + " is now fucking root");
                            //Step 8
                            for(Node n : sons)
                                send(n, new Message(frag, "NEW"));
                        }
                    }
                }

                //NEW
                else if(flag == "NEW") {
                    //Step 10
                    father = sender;
                    frag = (Integer) content;

                    //If it turns out another fragment was the root,
                    //remove it from sons
                    int pos = sons.indexOf(father);
                    if(pos != -1)
                        sons.remove(pos);
                    
                    for(Node s : sons)
                        send(s, new Message(content, flag));

                    if(sons.isEmpty()) {
                        send(father, new Message(null, "ECHO"));
                        phase++;
                    }
                    status = Status.ROOT;
                }

                //ECHO
                else if(flag == "ECHO") {
                    echoCounter++;
                    if(echoCounter == sons.size()) {
                        status = Status.BEGIN;
                        if(father != this) {
                            send(father, new Message(content, flag));
                            phase++;
                        }
                        else { //I am root
                            //Step 1
                            phase++;
                            processPulse();
                        }
                    }
                }

                //DISP
                else if(flag == "DISP") {
                    disp();
                }

                it.remove();
            }
            else {
                System.out.println("FURUETEU : " + getID());
                System.out.println(phase + " < " + mp.phase);
                System.out.println(m.getFlag());
            }
        }
    }

    /*
    @Override
    public void onMessage(Message message) {
    }*/

    @Override
    public void onSelection() {
        // what to do when this node is selected by the user
    }
}
