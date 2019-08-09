'''
    This script is created by Csenger Kovácsházi as a Thesis work in tha lab:
            Hungaryan Academy of Science, Research Centre of Natural Sciences
            Institute of Enzimology, Genome Stability Research Group
            Group leader: Dávid Szüts, PhD
        
    Function: Analyze Cas9 cutting pattern from NGS amplicon data
    Input: Needle file with the global alignment of the NGS reads one by one with the reference seqene
    Output: The alignments with the data of relevant events at the cutting position
    
    How it works:
    Read in the alignments from the input files
    With the alignment info, check if there is mismatch or gap at the cutting position
    If yes decide what is that: insertion, deletion, substitution or deletion with insertion
    Store relevant data and print out to a file.
'''

import re

tholerance = 2 # Says how far the indel can be from the cutting position
               # 2 is recommended
mismatch_minscore = 0 #What is the minimal score for mismatch what it must have to take as a mismatch (0 or higher)
substitution_minlength = 2 #What is the minimal length of substitution. (1 or higher)

class read: #Classs for reads we read from the file
    __slots__ = ['ID','name','ref','seq','aligninfo','score','has']
    #name - name of the read from the file
    #ref - reference sequence of the read (this can contain insertions)
    #seq - sequence of the read (this can contain deletions)
    #aligninfo - this contains if a base is a match, mismatch or gap
    #score - Score value of the alignment
    #has - The ID of the event what it has. First part. ID of the read, second part can be:
        #can be: 
        # - False -> doesn't have relevant event
        # - Error -> Relevand event but can't decide what is this
        # - D = deletion -> relevant deletion
        # - I = insertion -> relevant insertion with clear ends
        # - Imm =  ins with MM -> relevant insertion with mismathces at the side
        # - InD = insertion and deletion-> relevant deletion with mismatches at the side
        # - S = substitution -> mismatches at the cutting position
        # - Snd = substitution with 
        #     nonrelevant del -> also means substitution, but found in an other way (cause the read has gap not at the cutting positon)
        
            
    def printall(self): #print all the data of the read
        print("ID: ",self.ID)
        print("read name: "+self.name);
        print("reference  sequence: ",self.ref);
        print("aligninfo:           ",self.aligninfo)
        print("read sequence:       ",self.seq);
        print("Event ID: ",str(self.has))
        print("Score: ",str(self.score));
        print("_______________________________")
        return
    
    def give_line(self):
        line = self.ID+'\t'
        line += self.name+'\t'
        line += self.ref+'\t'
        line += self.aligninfo+'\t'
        line += self.seq+'\t'
        line += str(self.has)+'\t'
        line += str(self.score)+'\t'
        line += '\n'
        return line;
        
    
class deletion: #Class for event - deletions
    __slots__ = ['ID','readname','read','ref','startpos','endpos','mh_seq']
    
    def __init__(self,ID,read,startpos,endpos):
        self.ID = ID
        self.readname = read.name #name of the read what contains the deletion
        self.read = read.seq      #seq of the read
        self.ref = read.ref       #reference seq of the read
        self.startpos = startpos  #deletion start position
        self.endpos = endpos      #deletion end position
        self.mh_seq = common_start(self.read[self.endpos:], self.ref[startpos:endpos])
                                  #sequence of the microhomology what the deletion has
        
    def printall(self): #Print out all the data what we have from the deletion
        print("ID: ",self.ID)
        print("readname: ",self.readname);
        print("read seq:      ",self.read);
        print("reference seq: ",self.ref);
        print("startpos: ",self.startpos);
        print("endpos: ",self.endpos);
        print("MH sequence:"+self.mh_seq)
        print("_______________________________")
        return;
    def give_line(self):
        line = self.ID+'\t'
        line += self.readname+'\t'
        line += self.read+'\t'
        line += self.ref+'\t'
        line += str(self.startpos)+'\t'
        line += str(self.endpos)+'\t'
        if self.mh_seq == '':
            line += 'NA \t'
        else:
            line += self.mh_seq+'\t'
        line += '\n'
        return line

class insertion:
    __slots__ = ['ID','readname','read','ref','startpos','endpos','ins_seq']

    def __init__(self,ID,read,startpos,endpos):
        self.ID = ID
        self.readname = read.name  #name of the read what contains the insertion
        self.read = read.seq       #seq of the read
        self.ref = read.ref        #reference seq of the read
        self.startpos = startpos   #insertion start position
        self.endpos = endpos       #insertion end position
        
    def printall(self): #Print out all the data what we have from the insertion
        print("ID: ",self.ID)
        print("readname: ",self.readname);
        print("read seq:      ",self.read)
        print("reference seq: ",self.ref)
        print("startpos: ",self.startpos)
        print("endpos: ",self.endpos);
        print("_______________________________")
        return
    
    def give_line(self):
        line = self.ID+'\t'
        line += self.readname+'\t'
        line += self.read+'\t'
        line += self.ref+'\t'
        line += str(self.startpos)+'\t'
        line += str(self.endpos)+'\t'
        line += '\n'
        return line

class in_del:# class for event deletion with insertion
    __slots__ = ['ID','readname','read','ref','startpos','endpos','del_seq','ins_seq']
    
    def __init__(self,ID,read,start,end,realstart,realend,what):
        self.ID = ID
        self.readname = read.name  #name of the read what contains the indel
        self.read = read.seq       #seq of the read
        self.ref = read.ref        #reference seq of the read
        self.startpos = realstart     #The start position of the event; start = start position of the gap
        self.endpos = realend         #The end position of the event; end = end position of the gap
        
        if what == 'deletion':
            self.del_seq = read.ref[realstart:realend] #sequence of the deletion
            self.ins_seq = read.seq[realstart:start]+read.seq[end:realend] #sequence of the insertion
                                            #we can have mismatches at both ends of the deletion
        
        elif what == 'insertion':
            self.ins_seq = read.seq[realstart:realend]
            self.del_seq = read.ref[realstart:start]+read.ref[end:realend]
        
        else:
            print('ERROR in in-del',ID,'wrong \'what\' has been given')
            
    def give_line(self):
        line = self.ID+'\t'
        line += self.readname+'\t'
        line += self.read+'\t'
        line += self.ref+'\t'
        line += str(self.startpos)+'\t'
        line += str(self.endpos)+'\t'
        line += self.del_seq+'\t'
        line += self.ins_seq+'\t'
        line += '\n'
        return line
    
class substitution: #Class for the event substitution
    
    __slots__ = ['ID','readname','read','ref','startpos','endpos']
    
    def __init__(self,ID,read,startpos,endpos):
        self.ID = ID
        self.readname = read.name  #name of the read what contains substitution
        self.read = read.seq       #seq of the read
        self.ref = read.ref        #reference seq of the read
        self.startpos = startpos   #start position of the substitution
        self.endpos = endpos       #end position of the substitution
    
    def give_line(self):
        line = self.ID+'\t'
        line += self.readname+'\t'
        line += self.read+'\t'
        line += self.ref+'\t'
        line += str(self.startpos)+'\t'
        line += str(self.endpos)+'\t'
        line += '\n'
        return line
        
def common_start(sa, sb): # Check common start of two reads. Good for analysis of microhomologies
    def _iter():
        for a, b in zip(sa.upper(), sb.upper()):
            if a == b:
                yield a
            else:
                return
    return ''.join(_iter())

def openfile(path,MinScore): #path - path to the file; MinScore - minimal score what a read needs to reach

    file = None;
    while True: #Check the source path
        try:
            file = open(path,'r');
            break;
        except FileNotFoundError:
            print("the given path for the input file is wrong.")
            path = input("Please add the right path: ") 
    
        
    file.readline(); #read the first ###### row
    line = "null";
    
    while (not line.startswith("#===")): #Go thorough file header
        line = file.readline();
    
    alignments = []; #in this we'll store the alignments. This will be the output
    while (line != ""): #read till the end of the file
        
        if(line.startswith("#===")): # Here we jump through non-relevant lines and get info from the relevant ones
            newRead = read();
            
            file.readline();file.readline();file.readline(); #jump
                #name of the read______________
            newRead.name = file.readline().split()[2];
            
            #bigger jump
            while(not line.startswith("# Score")):
                line = file.readline();
            
                #score of the read_________________
            newRead.score = float(line.split()[-1]);
            
            if newRead.score < MinScore: #Jump to the next read if true
                del newRead;
                while(not line.startswith("#===")):
                      line = file.readline();
                while(not line.startswith("#===")):
                      line = file.readline();      
            else:
                newRead.ID = 'R'+str(len(alignments))
                #bigger jump
                while(not line.startswith("#===")):
                    line = file.readline();
                file.readline();
                #seq of the reference______________
            
                newRead.ref = file.readline();
                newRead.ref = newRead.ref.split()[2];
            
                #info about alignment___________
            
                newRead.aligninfo = file.readline()[21:-1]
                newRead.aligninfo = re.sub(" ", "_", newRead.aligninfo) #replace 'space' with _
            
                #seq of the read_____________
                newRead.seq = file.readline().split()[2];
            
                alignments.append(newRead); #store the read
                del newRead #memory optimalization
            
        line = file.readline();
        
    del line
    file.close()
    
    
    return alignments;

def mismatch(From,To,aligninfo):#Analyze of mismatches at the sides of the deletions.
                                #From/To: Determines the region and the direction where we look for mismatch.  
                                    #IF: From < To we're at the right side of the deletion -> go right
                                    #IF: From > To we're at the left side of the deletions, so go reverse direction
    #How it works:
        #Scores the alignment and take the start with the highest score
        #Scoring: +1 for mismatch, -1 for match, 0 for gap
    
    if From < To: #We go right
        scores = []
        
        for i in range(From,To): #Make scoring
            score = aligninfo[From-1:i].count('|')*-1 + aligninfo[From:i].count('.')*1
            scores.append(score)
        if scores == []: #From == To
            return From
        bestScore = max([(value,index) for index,value in enumerate(scores)]) #The value and index of the last max value in the list
        if bestScore[0] < mismatch_minscore:
            return From
        else:
            return From+bestScore[1]
        
    else: # From > To: #We go left
        invertinfo = aligninfo[To:From][::-1]                                   #Here happens the same like before
        scores = []                                                             #But now we wanna go backward, so we invert our seq

        for i in range(len(invertinfo)):                                        #And beacuse we have to invert we also strip it only to it's relevant part
            score = invertinfo[:i].count('|')*-1 + invertinfo[:i].count('.')*1
            scores.append(score)
        if scores == []:
            return From
        bestScore = max([(value,index) for index,value in enumerate(scores)]) #The value and index of the last max value in the list
        if bestScore[0] < mismatch_minscore:
            return From
        else:
            return From-bestScore[1]#Here we have a reverted scoretable because of the reverted aligninfo. That's why we need minus instead of plus

def check_substitution(read,cutpos): #Check if there is substitution around the cutting position (cutpos)
    rel_part = read.aligninfo[cutpos-tholerance:cutpos+tholerance] #we look substitution around the cutpos in range tholerance
    if rel_part.count('.') > 0: #If there is at least two mismatches, Check where the mismatching starts and ends
        position = cutpos-tholerance+rel_part.find('.')
        realstart = mismatch(position,0,read.aligninfo)
        realend = mismatch(position,len(read.aligninfo),read.aligninfo)
        if realend-realstart >= substitution_minlength:
            return [realstart,realend]
        else:
            return False
    else:
        return False

def findevent(aligns,cutting_position):
    dels = []
    ins = []
    mixed = []
    subst = []
    
    #Here we say that there is only one event what is at CP. When we have that, we check what kind of event is that
                #Than break the cycle, save the read in at the right place
                #If the event is not at cp -> check the following event
    
    for read in aligns: #iterate through all the reads
        events = list(re.finditer('_+',read.aligninfo)) #look for gap in aligninfo
        
        if not events: #If there is no gap, check and decide if it has substitution
            has_sub = check_substitution(read,cutting_position)
            if not has_sub:
                read.has = False
            else:
                start = has_sub[0]
                end = has_sub[1]
                ID = read.ID+'_S'+str(len(subst))
                read.has = ID
                subst.append(substitution(ID,read,start,end))        
        
        else: #If there is at least one gap
            placed = False
            for match in events:
                start = match.span()[0]
                end = match.span()[1]
                
                realstart = mismatch(start,0, read.aligninfo) # Is there mismatch before the gap?
                realend = mismatch(end,len(read.aligninfo),read.aligninfo) #Is there mismatch after the gap?
                
                
                if realstart <= cutting_position+tholerance and realend >= cutting_position-tholerance: #If this gap is relevant:     
                    #Check what is this                                                                   #Is around cutting position                   
                    if read.seq[match.span()[0]] == '-':
                        #The read seq has gap, so it is deletion based
                        if realstart == match.span()[0] and realend == match.span()[1]:
                            #Doesn't have mismatch -> clear deletion
                            ID = read.ID+'_D'+str(len(dels))
                            read.has = ID
                            dels.append(deletion(ID,read,start,end))
                            placed = True
                            break;
                        else:
                            #Has mismatch -> deletion with insertion
                            ID = read.ID+'_DwI'+str(len(mixed))
                            read.has = ID
                            mixed.append(in_del(ID,read,start,end,realstart,realend,'deletion'))
                            placed = True
                            break;
                      
                    elif read.ref[match.span()[0]] == '-':
                        #reference seq has gap, so it is insertion based
                        if realstart == match.span()[0] and realend == match.span()[1]:
                            #doesn't have mismatch -> clear insertion
                            ID = read.ID+'_I'+str(len(ins))
                            read.has = ID
                            ins.append(insertion(ID,read,realstart,realend))
                            placed = True
                            break;
                        else:
                            #has mismatch -> insertion with deletion
                            ID = read.ID+'_IwD'+str(len(mixed))
                            read.has = ID
                            mixed.append(in_del(ID,read,start,end,realstart,realend,'insertion'))
                            placed = True
                            break;
                              
                    else:
                        #something went wrong
                        ID = 'Error'
                        read.has = ID
                        print('NOTE: System can\'t analyze the read:',read.name)
                        placed = True
                        break;
                        
            if not placed: #Has gap but not relevant
                
                has_sub = check_substitution(read,cutting_position)
                if not has_sub: #Doesn't have substitution, so doesn't have anything
                    read.has = False
                else: #Has substitution
                    ID = read.ID+'_Snd'+str(len(subst))
                    read.has = ID
                    subst.append(substitution(ID,read,has_sub[0],has_sub[1]))
                    
                
    all_events = [dels,ins,mixed,subst]         
    return all_events


def deep_result_analysis(reads):
    noev = 0
    strange = 0
    dels = 0
    ins = 0
    mixes = 0
    subst = 0
    
    for i in reads:
        if i.has == False:
            noev += 1
        elif i.has == 'strange':
            strange += 1
        elif i.has == 'deletion':
            dels += 1
        elif i.has == 'insertion' or i.has == 'ins with MM':
            ins += 1
        elif i.has == 'del with MM':
            mixes += 1
        elif i.has == 'substitution' or i.has == 'subst with nonrel del':
            subst += 1
    print('noevent:',noev)
    print('strange:',strange)
    print('deletions:',dels)
    print('insertions:',ins)
    print('insertion and deletion:',mixes)
    print('substitution:',subst)
    
def initiate(path,cutting_position,minimal_score = 200):
    print('reading in file',path)
    reads = openfile(path,minimal_score) #returns a list with the reads
    events = findevent(reads,cutting_position) #returns a list of lists what contain all the types of events: 
    return [reads]+events                               #[deletions,insetrions,InsAndDels,Substitution]

def saveall(dataset,outfile = 'CAout.txt'):
    #save reads:_______________________________
    with open('Reads_'+outfile,'w') as O:
        O.write("ID \t Name \t Read_sequence \t Aligninfo \t  Reference_sequenc \t Has_event \t Alignment_Score \n")
        for read in dataset[0]:
            O.write(read.give_line())
        print('file with the reads has been created')
    
    with open('Dels_'+outfile,'w') as O:
        O.write("ID \t Read_name \t Read_sequence \t Reference_sequence \t Startp_position \t End_position \t Microhomology_sequence \n")
        for Del in dataset[1]:
            O.write(Del.give_line())
        print('file with the deletions has been created')
    
    with open('Ins_'+outfile,'w') as O:
        O.write("ID \t Read_name \t Read_sequence \t Reference_sequence \t Startp_position \t End_position \n")
        for ins in dataset[2]:
            O.write(ins.give_line())
        print('file with the insertions has been created')
    
    with open('Indels_'+outfile,'w') as O:
        O.write("ID \t Read_name \t Read_sequence \t Reference_sequence \t Startp_position \t End_position \t Deleted_sequence \t Inserted_sequence \n")
        for mix in dataset[3]:
            O.write(mix.give_line())
        print('file with the deletions and insertions has been created')
    
    with open('Subs_'+outfile,'w') as O:
        O.write("ID \t Read_name \t Read_sequence \t Reference_sequence \t Startp_position \t End_position \n")
        for sub in dataset[4]:
            O.write(sub.give_line())
        print('file with the substitutions has been created')
