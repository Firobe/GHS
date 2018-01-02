import jbotsim.*;
import jbotsim.Message;
import java.util.*;

public class GHSNode extends Node{
    public enum Status {
        INIT,
        BEGIN,
        READY,
        ROOT;
    }

    public class Invitation {
        public Node sender;
        public Integer frag;
        public Invitation(Node s, Integer f) {
            sender = s; frag = f;
        }
    }

    private Integer frag; //ID of my fragment
    private Node father; //Father in MST
    private ArrayList<Node> sons; //Sons in MST
    private Status status; //Current status (INIT, ..., ROOT)
    private Integer phase; 
    //Increased each time the root is updated during a merge
    private List<Node> neighbors;
    //List of neighbors in the graph
    //(to prevent multiple calls to getNeighbors)

    //Number of RFRAG, MCOE and ECHO to receive before propagating
    private Integer fragCounter;
    private Integer mcoeCounter;
    private Integer echoCounter;

    private Link curMCOE; //Current computed MCOE
    private Node goodboy; //Son that sent MCOE (null if irrelevant)

    private ArrayList<Invitation> mergeList;
    //List of invitations (MERGE requests) before READY status

    // Initialize the node variables
    @Override
    public void onStart() {
        frag = getID();
        father = this;
        sons = new ArrayList<Node>();
        mergeList = new ArrayList<Invitation>();
        status = Status.INIT;
        phase = 0;
        neighbors = getNeighbors();
    }

    // Called at the beginning of a new phase
    public void init() {
        mcoeCounter = 0;
        curMCOE = null;
        goodboy = null;
        echoCounter = 0;
        fragCounter = 0;
    }

    /**
     * Strict order between edges
     * (cf. poly)
     */
    private boolean greaterThan(Link l) {
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

    /**
     * Used at the end to display the computed MST in dot
     */
    private void disp() {
        for(Node s : sons) {
            System.out.println(getID() + "--"
                    + s.getID() + ";");
            send(s, new Message(null, "DISP"));
        }
    }

    /**
     * Called when the MCOE is completely computed
     * It must be either forwarded to father or (in case of the root)
     * trigger a ACK broadcast.
     *
     * If there was no MCOE, the algorithm is finished and we
     * display the MST.
     */
    private void sendMCOE() {
        if(father != this)
            send(father, new Message(curMCOE, "MCOE"));
        else { //I am root
            if(curMCOE == null) {//Step 4
                System.out.println("Algorithme terminé");
                System.out.println("graph {");
                disp();
            }
            else {
                //Step 5
                ackProcess(curMCOE);
                //Simulate sending a ACK with curMCOE to self
            }
        }
    }

    /**
     * Wrapper around a message to include
     * the phase of its sender when it was sent
     */
    public class MessNPhase {
        public Object object;
        public Integer phase;
        public MessNPhase(Object o, Integer p) {
            object = o;
            phase = p;
        }
    }

    /**
     * For each message sent, wrap our pahse number inside
     */
    @Override
    public void send(Node d, Message m) {
        MessNPhase mp = new MessNPhase(m.getContent(), phase);
        super.send(d, new Message(mp, m.getFlag()));
        System.out.println("sending " + m.getFlag() + " from " + getID() + " to " + d.getID());
    }

    /**
     * Called when a PULSE is received or to begin a new phase
     */
    private void processPulse() {
        //Step 1
        init();
        for(Node n : neighbors) //Send PULSE to sons and FRAG to others
            if(sons.contains(n))
                send(n, new Message(null, "PULSE"));
            else if (n != father){
                send(n, new Message(null, "FRAG"));
                fragCounter++; //Remember we must receive a RFRAG
            }
    }

    /**
     * Called when ACK is received
     * Re-orientate the tree and send MERGE requests
     * if necessary
     */
    private void ackProcess(Link rMCOE) {
        //Step 5

        //Propagate ACK to sons
        for(Node n : sons)
            send(n, new Message(rMCOE, "ACK"));

        //Update frag number
        frag = rMCOE.endpoint(0).getID();
        status = Status.READY;

        //If on the path to MCOE, change orientation
        if(curMCOE != null && rMCOE.equals(curMCOE)) { 
            if(father != this) //Not current root
                sons.add(father);
            father = goodboy;
            if(goodboy != null) //If not new root, swap
                sons.remove(sons.indexOf(goodboy));
            else { //If I am the one who computed the MCOE, MERGE
                father = this;
                //Step 6
                send(curMCOE.endpoint(1), new Message(
                            frag, "MERGE"));
            }
        } 

        //Step 7
        //Process stored MERGE requests
        for(Invitation p : mergeList) {
            sons.add(p.sender);
            //Systematically add to sons, although may be removed later
            if(father == this //If I'm the fragment root
                && rMCOE.equals(new Link(this, p.sender))
                //And the edge between me and the sender is my MCOE
                && frag < p.frag){ //And my fragment is inferior
                    status = Status.ROOT; //Then I'm the new root
                    phase++;
                }
        }

        //If I'm the root, warn my sons with NEW
        if(status == Status.ROOT) {
            System.err.println(getID() + " is now a root");
            //Step 8
            for(Node n : sons)
                send(n, new Message(frag, "NEW"));
        }
        //Update curMCOE in case someone wants to merge with me later
        curMCOE = rMCOE; 
    }

    @Override
    public void onClock() {
        List<Message> received = getMailbox();
        //Initial pulse for the first phase
        if(status == Status.INIT) {
            processPulse();
            status = Status.BEGIN;
        }

        //Reading messages
        for(Iterator<Message> it = received.iterator() ; it.hasNext() ; ) {
            //Unwrap the sender's phase from the message
            Message m = it.next();
            MessNPhase mp = (MessNPhase) m.getContent();
            //Only process a message if its not from the future
            if(mp.phase <= phase || m.getFlag() != "FRAG") {
                //Create aliases
                String flag = m.getFlag();
                Object content = mp.object;
                Node sender = m.getSender();

                //Now filter the name of the message and proceed
                //accordingly

                //PULSE
                if(flag == "PULSE"){
                    phase = mp.phase;
                    processPulse();
                }
	   
                //FRAG
                else if(flag == "FRAG") //Step 2
                    send(sender, new Message(frag, "RFRAG"));
                
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
                    Link l = (Link) content;
                    //Inside step 3
                    if(sender != null) //If MCOE is not from myself
                        mcoeCounter++;

                    //Update local minimum
                    if(sender != null && greaterThan(l)) {
                        curMCOE = l;
                        goodboy = sender; //Remember who sent it
                    }

                    //Received MCOE from every son, send to father
                    if(mcoeCounter == sons.size() && fragCounter == 0)
                        sendMCOE();
                }

                //ACK
                else if(flag == "ACK") {
                    ackProcess((Link) content);
                }

                //MERGE
                else if(flag == "MERGE") {
                    if(status == Status.BEGIN
                            && phase == mp.phase) {
                        //Step 7
                        //If I'm not READY, store the request for later
                        mergeList.add(new Invitation(sender,
                                    (Integer) content));
                    }
                    else if(status == Status.ROOT
                            || phase > mp.phase) {//
                        //Step 9
                        //If the root has already been decided
                        //Warns the late merge-wannabe with NEW
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
                            phase++;
                        }
                        if(status == Status.ROOT) {
                            System.out.println(getID() + " is now a root");
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
		    phase = mp.phase;
                    //If it turns out another fragment was the root,
                    //remove it from sons
                    int pos = sons.indexOf(father);
                    if(pos != -1)
                        sons.remove(pos);
                    
                    for(Node s : sons)
                        send(s, new Message(content, flag));

                    if(sons.isEmpty()) {
                        status = Status.BEGIN;
                        mergeList.clear();
                        send(father, new Message(null, "ECHO"));
                    }
                }

                //ECHO
                else if(flag == "ECHO") {
                    if(status != Status.BEGIN) {
                        echoCounter++;
                        if(echoCounter == sons.size()) {
                            status = Status.BEGIN;
                            mergeList.clear();
                            if(father != this) {
                                send(father, new Message(content, flag));
                }
                            else { //I am root
                                //Step 1
                                processPulse();
                            }
                        }
                    }
                    else {
                        send(sender, new Message(null, "PULSE"));
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
