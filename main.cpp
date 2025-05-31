#include <bits/stdc++.h>
#include <fstream>

using namespace std;

class Node;

class Element;

class Node {
private:
    string name;
    double voltage;
    vector<Element *> connected_elements;
public:
    Node(string name_) : name(name_) {}

    vector<Element *> getConnectedElements() { return connected_elements; }

    string getName() { return name; }

    void addConnectedElement(Element *e) {
        connected_elements.push_back(e);
    }

    void removeConnectedElement(Element *e) {
        auto x = find(connected_elements.begin(), connected_elements.end(), e);
        if (x != connected_elements.end())
            connected_elements.erase(x);
    }

    void setName(string new_name) {
        name = new_name;
    }

    double getVoltage() { return voltage; }
};

class Circuit {
protected:
    vector<Element *> elements;
    vector<Node *> nodes;
    map<string, Node *> node_access;
    Node *gnd=nullptr;
public:
    vector<Element *> &getElements() { return elements; }

    map<string, Node *> getNodeAccess() { return node_access; }

    Node *getCreateNode(string name) {
        if (!node_access.count(name)) {
            Node *newnode = new Node(name);
            nodes.push_back(newnode);
            node_access[name] = newnode;
        }
        return node_access[name];
    }

    void addElement(Element *element) {
        elements.push_back(element);
    }

    Node* getGround() {return gnd;}
    vector<Node *> getNodes() {return nodes;}

    bool setGround(string name){
        if(!node_access.count(name))
            return false;
        gnd=node_access[name];
        return true;
    }

    bool removeGround(string name){
        if(gnd!=node_access[name])
            return false;
        gnd=nullptr;
        return true;
    }

    bool renameNode(string old_name, string new_name, string& err) {
        if (node_access.find(old_name) == node_access.end()) {
            err = "ERROR: Node " + old_name + " does not exist in the circuit\n";
            return false;
        }
        if (node_access.find(new_name) != node_access.end()) {
            err = "ERROR: Node name " + new_name + " already exists\n";
            return false;
        }
        Node* node = node_access[old_name];
        node->setName(new_name);
        node_access.erase(old_name);
        node_access[new_name] = node;
        err = "SUCCESS: Node renamed from " + old_name + " to " + new_name + "\n";
        return true;
    }

};

class Element {
protected:
    string name;
    double value;
    Node *node1;
    Node *node2;
public:
    Element(string name_, double value_, Node *n1, Node *n2) :
            name(name_), value(value_), node1(n1), node2(n2) {}

    string getName() {
        return name;
    }

    virtual string getType() = 0;

    double getValue() { return value; }

    Node *getFirstNode() { return node1; }

    Node *getSecondNode() { return node2; }
    virtual double getCurrent() {

    }

    virtual double getVoltage() {

    }

    virtual ~Element() = default;
};

class Resistor : public Element {
public:
    Resistor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Resistor"; }

    double getCurrent() override {
    }

    double getVoltage() override {

    }
};

class Capacitor : public Element {
public:
    Capacitor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Capacitor"; }

    double getCurrent() override {

    }

    double getVoltage() override {

    }
};

class Inductor : public Element {
public:
    Inductor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Inductor"; }

    double getCurrent() override {

    }

    double getVoltage() override {

    }
};

class Diode : public Element {
protected:
    string model;
public:
    Diode(string name, Node *n1, Node *n2) :
            Element(name, 0, n1, n2), model("D") {}

    string getType() override { return "Diode"; }

    virtual double calculateCurrent(double voltage) {

    }
};

class ZenerDiode : public Diode {
private:
    double breakdownvoltage;
public:
    ZenerDiode(string name, Node *n1, Node *n2, double breakdownvoltage_ = 5)
            : Diode(name, n1, n2), breakdownvoltage(breakdownvoltage_) {
        model = "Z";
    }

    string getType() override { return "Zener Diode"; }

    double calculateCurrent(double voltage) override {

    }
};

class VoltageSource : public Element {
public:
    VoltageSource(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    virtual double getVoltage(double time) = 0;

    string getType() override { return "Voltage Source"; }
};

class DCVoltageSource : public VoltageSource {
public:
    DCVoltageSource(string name_, double value_, Node *n1, Node *n2)
            : VoltageSource(name_, value_, n1, n2) {}

    double getVoltage(double time) override {
        return getValue();
    }

    string getType() override { return "DC Voltage Source"; }
};

class SinusoidalVoltageSource : public VoltageSource {
private:
    double frequency;
    double offset;

public:
    SinusoidalVoltageSource(string name_, double amplitude, double frequency_, double offset_, Node *n1, Node *n2)
            : VoltageSource(name_, amplitude, n1, n2), frequency(frequency_), offset(offset_) {}

    double getVoltage(double time) override {
        return getValue() * sin(2 * M_PI * frequency * time) + offset;
    }

    double getFrequency(){return frequency;}

    double getOffset(){return offset;}

    string getType() override { return "Sinusoidal Voltage Source"; }
};

class PulseVoltageSource : public VoltageSource {
private:
    double V1;
    double V2;
    double TD;
    double TR;
    double TF;
    double TOn;
    double period;

public:
    PulseVoltageSource(string name_, double V1_, double V2_, double TD_, double TR_, double TF_, double TOn_, double period_, Node *n1, Node *n2)
            : VoltageSource(name_, V1_, n1, n2), V1(V1_), V2(V2_), TD(TD_), TR(TR_), TF(TF_), TOn(TOn_), period(period_) {}

    double getVoltage(double time) override {
        double t = fmod(time, period);

        if (t < TD) {
            return V1;
        }
        else if (t < TD + TR) {
            return V1 + (V2 - V1) * (t - TD) / TR;
        }
        else if (t < TD + TR + TOn) {
            return V2;
        }
        else if (t < TD + TR + TOn + TF) {
            return V2 - (V2 - V1) * (t - TD - TOn - TR) / TF;
        }
        else {
            return V1;
        }
    }

    string getType() override { return "Pulse Voltage Source"; }
};

class CurrentSource : public Element {
public:
    CurrentSource(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    virtual double getCurrent(double time) = 0;

    string getType() override { return "Current Source"; }
};

class DCCurrentSource : public CurrentSource {
public:
    DCCurrentSource(string name_, double value_, Node *n1, Node *n2)
            : CurrentSource(name_, value_, n1, n2) {}

    double getCurrent(double time) override {
        return getValue();
    }

    string getType() override { return "DC Current Source"; }
};

class SinusoidalCurrentSource : public CurrentSource {
private:
    double frequency;
    double offset;

public:
    SinusoidalCurrentSource(string name_, double amplitude, double frequency_, double offset_, Node *n1, Node *n2)
            : CurrentSource(name_, amplitude, n1, n2), frequency(frequency_), offset(offset_) {}

    double getCurrent(double time) override {
        return getValue() * sin(2 * M_PI * frequency * time) + offset;
    }

    double getFrequency(){return frequency;}

    double getOffset(){return offset;}

    string getType() override { return "Sinusoidal Current Source"; }
};

class PulseCurrentSource : public CurrentSource {
private:
    double I1;
    double I2;
    double TD;
    double TR;
    double TF;
    double TOn;
    double period;

public:
    PulseCurrentSource(string name_, double I1_, double I2_, double TD_, double TR_, double TF_, double TOn_, double period_, Node *n1, Node *n2)
            : CurrentSource(name_, I1_, n1, n2), I1(I1_), I2(I2_), TD(TD_), TR(TR_), TF(TF_), TOn(TOn_), period(period_) {}

    double getCurrent(double time) override {
        double t = fmod(time, period);

        if (t < TD) {
            return I1;
        }
        else if (t < TD + TR) {
            return I1 + (I2 - I1) * (t - TD) / TR;
        }
        else if (t < TD + TR + TOn) {
            return I2;
        }
        else if (t < TD + TR + TOn + TF) {
            return I2 - (I2 - I1) * (t - TD - TOn - TR) / TF;
        }
        else {
            return I1;
        }
    }

    string getType() override { return "Pulse Current Source"; }
};

class VCVS : public VoltageSource {
private:
    double gain;
    Node *controlNode1;
    Node *controlNode2;
public:
    VCVS(string name_, double gain_, Node *n1, Node *n2, Node *controlNode1_, Node *controlNode2_) :
            VoltageSource(name_, 0, n1, n2), gain(gain_), controlNode1(controlNode1_), controlNode2(controlNode2_) {}

    Node* getControlNode1() { return controlNode1; }
    Node* getControlNode2() { return controlNode2; }
    double getGain() { return gain; }

    double getVoltage(double time) override {
        double controlVoltage = controlNode1->getVoltage() - controlNode2->getVoltage();
        return gain * controlVoltage;
    }

    string getType() override {
        return "VCVS";
    }
};

class CCVS : public VoltageSource {
private:
    double gain;
    Element *controlElement;
public:
    CCVS(string name_, double gain_, Node *n1, Node *n2, Element *controlElement_)
            : VoltageSource(name_, 0, n1, n2), gain(gain_), controlElement(controlElement_) {}
    double getGain() { return gain; }
    double getVoltage(double time) override {
        double controlCurrent = controlElement->getCurrent();
        return gain * controlCurrent;
    }

    Element *getcontrolElement() {return controlElement;}

    string getType() override {
        return "CCVS";
    }
};

class VCCS : public CurrentSource {
private:
    double gain;
    Node *controlNode1;
    Node *controlNode2;
public:
    VCCS(string name_, double gain_, Node *n1, Node *n2, Node *controlNode1_, Node *controlNode2_) :
            CurrentSource(name_, 0, n1, n2), gain(gain_), controlNode1(controlNode1_), controlNode2(controlNode2_) {}

    Node* getControlNode1() { return controlNode1; }
    Node* getControlNode2() { return controlNode2; }
    double getGain() { return gain; }

    double getCurrent(double time) override {
        double controlVoltage = controlNode1->getVoltage() - controlNode2->getVoltage();
        return gain * controlVoltage;
    }

    string getType() override {
        return "VCCS";
    }
};

class CCCS : public CurrentSource {
private:
    double gain;
    Element *controlElement;
public:
    CCCS(string name_, double gain_, Node *n1, Node *n2, Element *controlElement_)
            : CurrentSource(name_, 0, n1, n2), gain(gain_), controlElement(controlElement_) {}

    double getGain() { return gain; }

    double getCurrent(double time) override {
        double controlCurrent = controlElement->getCurrent();
        return gain * controlCurrent;
    }

    Element *getcontrolElement() {return controlElement;}

    string getType() override {
        return "CCCS";
    }
};

class Controller {
public:
    string handleError(string name, Circuit *circuit, double value) {
        string type;
        switch (name[0]) {
            case 'R':
                type = "Resistor";
                break;
            case 'C':
                if(name.size()==0)
                    type = "Capacitor";
                break;
            case 'L':
                type = "Inductor";
                break;
            case 'D':
                type = "Diode";
                break;
        }
        if (value <= 0 && (type == "Resistor" || type == "Capacitor" || type == "Inductor"))
            return "Error: " + type + " cannot be zero or negative\n";
        for (auto element: circuit->getElements()) {
            if (!element)
                continue;
            if (element->getName() == name) {
                switch (name[0]) {
                    case 'R':
                        return "Resistor " + name + " already exists in the circuit\n";
                    case 'C':
                        return "Capacitor " + name + " already exists in the circuit\n";
                    case 'L':
                        return "Inductor " + name + " already exists in the circuit\n";
                    case 'D':
                        return "Diode" + name + " already exists in the circuit\n";
                    default:
                        return "Source already exists in the circuit\n";
                }
            }
        }
        return "";
    }

    string addNewElement(string node1, string node2, string name, double value, Circuit *circuit) {
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        switch (name[0]) {
            case 'R':
                element = new Resistor(name, value, n1, n2);
                break;
            case 'L':
                element = new Inductor(name, value, n1, n2);
                break;
            case 'C':
                element = new Capacitor(name, value, n1, n2);
                break;
        }
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string removeElement(string name, Circuit *circuit) {
        auto &elements = circuit->getElements();
        for (auto it = elements.begin(); it != elements.end(); it++) {
            if ((*it)->getName() == name) {
                (*it)->getFirstNode()->removeConnectedElement(*it);
                (*it)->getSecondNode()->removeConnectedElement(*it);
                Element *deleting = *it;
                elements.erase(it);
                delete deleting;
                return name + " removed successfully!\n";
            }
        }
        switch (name[0]) {
            case 'R':
                return "Error: Cannot delete resistor; component not found\n";
            case 'L':
                return "Error: Cannot delete inductor; component not found\n";
            case 'C':
                return "Error: Cannot delete capacitor; component not found\n";
            case 'D':
                return "Error: Cannot delete diode; component not found\n";
        }
    }

    string addDiode(string node1, string node2, string name, string model, Circuit *circuit) {
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        if (model == "Z") {
            ZenerDiode *zenerDiode = new ZenerDiode(name, n1, n2);
            circuit->addElement(zenerDiode);
            n1->addConnectedElement(zenerDiode);
            n2->addConnectedElement(zenerDiode);
            return "Zener Diode added successfully!\n";
        } else if (model == "D") {
            Diode *diode = new Diode(name, n1, n2);
            circuit->addElement(diode);
            n1->addConnectedElement(diode);
            n2->addConnectedElement(diode);
            return "Diode added successfully!\n";
        }
    }

    string addGround(string node, Circuit* circuit){
        if(circuit->getGround())
            return "Error! ground node already exists in the circuit\n";
        Node *n1 = circuit->getCreateNode(node);
        circuit->setGround(node);
        return "Ground added to node " + node +" successfully!\n";
    }

    string deleteGround(string node, Circuit* circuit){
        if(!circuit->getNodeAccess().count(node))
            return "Node does not exist\n";
        if(!circuit->removeGround(node))
            return "Error! Node is not set as ground\n";
        return "Ground removed from node "+node+" successfully!\n";
    }

    void showNodes(Circuit* circuit){
        cout<<"Available nodes:\n";
        if(!circuit->getNodeAccess().size())
            cout<<"No node exists in the circuit\n";
        int index=1;
        for(auto node: circuit->getNodeAccess()){
            cout<<node.first;
            if(index!=circuit->getNodeAccess().size())
                cout<<", ";
            else
                cout<<"\n";
            index++;
        }
    }

    void list(Circuit* circuit){
        cout<<"Elements:\n";
        if(!circuit->getElements().size())
            cout<<"No elements exists in the circuit\n";
        for(auto element: circuit->getElements()){
            cout<<element->getType()<<": "<<element->getName()<<", value: "<<element->getValue()<<endl;
        }
    }

    void listElement(string type, Circuit* circuit){
        switch(type[0]){
            case 'R':
                type = "Resistor";
                break;
            case 'C':
                type = "Capacitor";
                break;
            case 'L':
                type = "Inductor";
                break;
            case 'D':
                type = "Diode";
                break;
        }
        cout<<type<<"(s):\n";
        int index=0;
        for(auto element: circuit->getElements()){
            if(element->getType()==type) {
                cout << element->getType() << ": " << element->getName() << ", value: " << element->getValue()<< endl;
                index++;
            }
        }if(!index)
            cout<<"No "<<type<<" exist in the circuit\n";
    }

    string addDCVoltageSource(string name, string node1, string node2, double value, Circuit *circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new DCVoltageSource(name, value, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addSinusoidalVoltageSource(string name, string node1, string node2, double amplitude,
                                      double frequency, double offset, Circuit *circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new SinusoidalVoltageSource(name, amplitude, frequency, offset, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addPulseVoltageSource(string name, string node1, string node2, double V1, double V2, double TD, double TR,
                                 double TF, double TOn, double period, Circuit* circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new PulseVoltageSource(name, V1, V2, TD, TR, TF, TOn, period, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addDCCurrentSource(string name, string node1, string node2, double value, Circuit *circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new DCCurrentSource(name, value, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addSinusoidalCurrentSource(string name, string node1, string node2, double amplitude,
                                      double frequency, double offset, Circuit *circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new SinusoidalCurrentSource(name, amplitude, frequency, offset, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addPulseCurrentSource(string name, string node1, string node2, double I1, double I2, double TD, double TR,
                                 double TF, double TOn, double period, Circuit* circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new PulseCurrentSource(name, I1, I2, TD, TR, TF, TOn, period, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addVCVS(string name, string node1, string node2, string controlNode1, string controlNode2, double gain, Circuit* circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Node *cn1 = circuit->getCreateNode(node1);
        Node *cn2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new VCVS(name, gain, n1, n2, cn1, cn2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addCCVS(string name, string node1, string node2, string controlElement, double gain, Circuit* circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        auto elements = circuit->getElements();
        Element *element = nullptr;
        Element *cElement = nullptr;
        for (auto& e : elements){
            if (e->getName() == controlElement){
                cElement = e;
                break;
            }
        }
        element = new CCVS(name, gain, n1, n2, cElement);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addVCCS(string name, string node1, string node2, string controlNode1, string conrolNode2, double gain, Circuit* circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Node *cn1 = circuit->getCreateNode(node1);
        Node *cn2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new VCCS(name, gain, n1, n2, cn1, cn2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addCCCS(string name, string node1, string node2, string controlElement, double gain, Circuit* circuit){
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        auto elements = circuit->getElements();
        Element *element = nullptr;
        Element *cElement = nullptr;
        for (auto& e : elements){
            if (e->getName() == controlElement){
                cElement = e;
                break;
            }
        }
        element = new CCCS(name, gain, n1, n2, cElement);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    void showCircuitDetails(Circuit *circuit) {
        cout << "Circuit Details:\nElements:\n";
        for (auto element: circuit->getElements()) {
            cout << "type: " << element->getType() << "   name: " << element->getName() << "  value: "
                 << element->getValue()
                 << "  node1: " << element->getFirstNode()->getName() << "  node2: "
                 << element->getSecondNode()->getName() << endl;
        }
        cout << "Node Details:\n";
        for (auto node: circuit->getNodeAccess()) {
            cout << node.first << " : " << endl;
            if (node.second->getConnectedElements().size() == 0)
                cout << "No element is connected to this node!\n";
            else {
                for (auto element: node.second->getConnectedElements())
                    cout << element->getName() << endl;
            }
        }
        cout<<"GND:\n";
        if(circuit->getGround())
            cout<<circuit->getGround()->getName()<<endl;
        else
            cout<<"No ground is specified for this circuit\n";
    }bool isConnected(Circuit* circuit) {
        if (circuit->getNodes().empty()) return true;

        unordered_map<string, vector<string>> adjacency;
        for (Element* e : circuit->getElements()) {
            if (!e) continue;
            string a = e->getFirstNode()->getName();
            string b = e->getSecondNode()->getName();
            adjacency[a].push_back(b);
            adjacency[b].push_back(a);
        }

        unordered_set<string> visited;
        queue<string> q;
        string start = circuit->getNodes()[0]->getName();
        q.push(start);
        visited.insert(start);

        while (!q.empty()) {
            string current = q.front(); q.pop();
            for (const string& neighbor : adjacency[current]) {
                if (!visited.count(neighbor)) {
                    visited.insert(neighbor);
                    q.push(neighbor);
                }
            }
        }

        return visited.size() == circuit->getNodeAccess().size();
    }bool preAnalysisErrs(Circuit *circuit){
        if(!circuit->getGround()){
            cout<<"Error: No ground node detected in the circuit.\n";
            return false;
        }if(!isConnected(circuit)) {
            cout << "Error: Not all nodes are connected!\n";
            return false;
        }
        for(auto element: circuit->getElements()) {
            if (element->getType() == "CCCS") {
                auto *cc = dynamic_cast<CCCS *>(element);
                if (cc) {
                    Element *ctrl = cc->getcontrolElement();
                    bool found = false;
                    for (auto e : circuit->getElements()) {
                        if (e == ctrl) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        cout << "Error: Dependent source " << cc->getName()
                             << "  has an undefined control element.\n";
                        return false;
                    }
                }
            }
            if (element->getType() == "CCVS") {
                auto *cv = dynamic_cast<CCVS *>(element);
                if (cv) {
                    Element *ctrl = cv->getcontrolElement();
                    bool found = false;
                    for (auto e : circuit->getElements()) {
                        if (e == ctrl) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        cout << "Error: Dependent source " << cv->getName()
                             << "  has an undefined control element.\n";
                        return false;
                    }
                }
            }if (element->getType() == "VCCS") {
                auto *vc = dynamic_cast<VCCS *>(element);
                if (vc) {
                    Node *ctrl = vc->getControlNode1();
                    Node*ctrl2 = vc->getControlNode2();
                    bool found = false;
                    bool found2 = false;
                    for (auto e : circuit->getNodes()) {
                        if (e == ctrl) {
                            found = true;
                        }
                        if(e==ctrl2)
                            found2=true;
                        if(found && found2)
                            break;
                    }
                    if (!found || !found2 || (ctrl->getConnectedElements().size()==0) || (ctrl2->getConnectedElements().size()==0)) {
                        cout << "Error: Dependent source " << vc->getName()
                             << "  has an undefined control element.\n";
                        return false;
                    }
                }
            }if (element->getType() == "VCVS") {
                auto *vv = dynamic_cast<VCVS *>(element);
                if (vv) {
                    Node *ctrl = vv->getControlNode1();
                    Node*ctrl2 = vv->getControlNode2();
                    bool found = false;
                    bool found2 = false;
                    for (auto e : circuit->getNodes()) {
                        if (e == ctrl) {
                            found = true;
                        }
                        if(e==ctrl2)
                            found2=true;
                        if(found && found2)
                            break;
                    }
                    if (!found || !found2 || (ctrl->getConnectedElements().size()==0) || (ctrl2->getConnectedElements().size()==0)) {
                        cout << "Error: Dependent source " << vv->getName()
                             << "  has an undefined control element.\n";
                        return false;
                    }
                }
            }
        }
        cout<<"Seems fine!\n";
        return true;
    }

};

class View {
private:
    Circuit *circuit;
    Controller controller;
public:
    void run() {
        circuit = new Circuit();
        string input;
        smatch match;
        regex add_element(
                R"(^\s*add\s+([A-Za-z])(\w+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s*$)");
        regex remove_element(R"(^\s*delete\s+([A-Za-z])(\w+)\s*$)");
        regex add_diode(R"(^\s*add\s+(D)(\w+)\s+(\S+)\s+(\S+)\s+(D|Z)\s*$)");
        regex show_details(R"(^\s*show\s+details\s*$)");
        regex add_ground(R"(^\s*add\s+(\w+)\s+(\w+)\s*$)");
        regex delete_ground(R"(^\s*delete\s+(\w+)\s+(\w+)\s*$)");
        regex show_nodes(R"(^\s*.nodes\s*$)");
        regex list(R"(^\s*.list\s*$)");
        regex list_element(R"(^\s*.list\s+(\w+)\s*$)");
        regex rename_node(R"(^\s*\.rename\s+node\s+(\S+)\s+(\S+)\s*$)");
        regex new_file(R"(NewFile ([a-zA-Z0-9\\-_\\.:\\/()\\s]+))");
        regex addDCSource(
                R"(add\s+(VoltageSource|CurrentSource)(\S+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s*$)");
        regex addSINSource(
                R"(^\s*add\s+([A-Za-z])(\w+)\s+(\S+)\s+(\S+)\s+SIN\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s*$)");
        regex addPULSESource(
                R"(^\s*add\s+([A-Za-z])(\w+)\s+(\S+)\s+(\S+)\s+PULSE\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s*$)");
        regex addVCVS(
                R"(^\s*add\s+E(\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?[GMKkmnp]?)\s*$)");
        regex addVCCS(
                R"(^\s*add\s+G(\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?[GMKkmnp]?)\s*$)");
        regex addCCVS(R"(^\s*add\s+H(\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?[GMKkmnp]?)\s*$)");
        regex addCCCS(R"(^\s*add\s+F(\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?[GMKkmnp]?)\s*$)");
        regex preanalysischeck(R"(^\s*pre-analysis\s+check\s*$)");
        regex exit(R"(^exit$)");
        while (true) {
            getline(cin, input);
            if (regex_match(input, match, addDCSource)) {
                string name = match[2].str();
                string node1 = match[3].str();
                string node2 = match[4].str();
                string number = match[5].str();
                string unit = match[6].str();
                double value;
                try {
                    value = stod(number);
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value\n";
                    continue;
                }
                string error = controller.handleError(name, circuit, value);
                if (error != "") {
                    cout << error;
                    continue;
                }
                if (!unit.empty()) {
                    switch (unit[0]) {
                        case 'G':
                            value *= 1e9;
                            break;
                        case 'M':
                            value *= 1e6;
                            break;
                        case 'k':
                        case 'K':
                            value *= 1e3;
                            break;
                        case 'm':
                            value *= 1e-3;
                            break;
                        case 'u':
                            value *= 1e-6;
                            break;
                        case 'n':
                            value *= 1e-9;
                            break;
                    }
                }
                if (match[1] == "VoltageSource")
                    cout << controller.addDCVoltageSource(name, node1, node2, value, circuit);
                else
                    cout << controller.addDCCurrentSource(name, node1, node2, value, circuit);
            } else if (regex_match(input, match, addSINSource)) {
                if (match[1] != "V" && match[1] != "I") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                string name = match[1].str() + match[2].str();
                string node1 = match[3].str();
                string node2 = match[4].str();
                string number[3] = {match[5].str(), match[7].str(), match[9].str()};
                string unit[3] = {match[6].str(), match[8].str(), match[10].str()};
                double value[3];
                int check = 0;
                for (int i = 0; i < 3; i++) {
                    try {
                        value[i] = stod(number[i]);
                    } catch (const invalid_argument &e) {
                        cout << value[i] << endl;
                        cout << "Error: Invalid numeric value\n";
                        check = 1;
                        break;
                    }
                }
                if (check)
                    continue;
                check = 0;
                for (int i = 0; i < 3; i++) {
                    string error = controller.handleError(name, circuit, value[i]);
                    if (error != "") {
                        cout << error;
                        check = 1;
                        break;
                    }
                }
                if (check)
                    continue;
                for (int i = 0; i < 3; i++) {
                    if (!unit[i].empty()) {
                        switch (unit[i][0]) {
                            case 'G':
                                value[i] *= 1e9;
                                break;
                            case 'M':
                                value[i] *= 1e6;
                                break;
                            case 'k':
                            case 'K':
                                value[i] *= 1e3;
                                break;
                            case 'm':
                                value[i] *= 1e-3;
                                break;
                            case 'u':
                                value[i] *= 1e-6;
                                break;
                            case 'n':
                                value[i] *= 1e-9;
                                break;
                        }
                    }
                }
                if (value[2] <= 0) {
                    cout << "Frequency can not be Zero or Negative!\n";
                    continue;
                }
                if (match[1] == "V")
                    cout << controller.addSinusoidalVoltageSource(name, node1, node2, value[1], value[2], value[0],
                                                                  circuit);
                else
                    cout << controller.addSinusoidalCurrentSource(name, node1, node2, value[1], value[2], value[0],
                                                                  circuit);
            } else if (regex_match(input, match, addPULSESource)) {
                if (match[1] != "V" && match[1] != "I") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                string name = match[1].str() + match[2].str();
                string node1 = match[3].str();
                string node2 = match[4].str();
                string number[7] = {match[5].str(), match[7].str(), match[9].str(), match[11].str(), match[13].str(),
                                    match[15].str(), match[17].str()};
                string unit[7] = {match[6].str(), match[8].str(), match[10].str(),
                                  match[12].str(), match[14].str(), match[16].str(), match[18].str()};
                double value[7];
                int check = 0;
                for (int i = 0; i < 7; i++) {
                    try {
                        value[i] = stod(number[i]);
                    } catch (const invalid_argument &e) {
                        cout << value[i] << endl;
                        cout << "Error: Invalid numeric value\n";
                        check = 1;
                        break;
                    }
                }
                if (check)
                    continue;
                check = 0;
                for (int i = 0; i < 7; i++) {
                    string error = controller.handleError(name, circuit, value[i]);
                    if (error != "") {
                        cout << error;
                        check = 1;
                        break;
                    }
                }
                if (check)
                    continue;
                for (int i = 0; i < 7; i++) {
                    if (!unit[i].empty()) {
                        switch (unit[i][0]) {
                            case 'G':
                                value[i] *= 1e9;
                                break;
                            case 'M':
                                value[i] *= 1e6;
                                break;
                            case 'k':
                            case 'K':
                                value[i] *= 1e3;
                                break;
                            case 'm':
                                value[i] *= 1e-3;
                                break;
                            case 'u':
                                value[i] *= 1e-6;
                                break;
                            case 'n':
                                value[i] *= 1e-9;
                                break;
                        }
                    }
                }
                if (value[2] < 0 || value[3] < 0 || value[4] < 0 || value[5] < 0 || value[6] < 0) {
                    cout << "Time parameters can not be Negative!\n";
                    continue;
                }
                if (value[6] < value[5] + value[4] + value[3] + value[2] + value[1]) {
                    cout << "Period shorter than enough!\n";
                    continue;
                }
                if (match[1] == "V")
                    cout << controller.addPulseVoltageSource(name, node1, node2, value[0], value[1], value[2], value[3],
                                                             value[4], value[5], value[6], circuit);
                else
                    cout << controller.addPulseCurrentSource(name, node1, node2, value[0], value[1], value[2], value[3],
                                                             value[4], value[5], value[6], circuit);
            } else if (regex_match(input, match, add_element)) {
                if (match[1] != "R" && match[1] != "L" && match[1] != "C") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                string name = match[1].str() + match[2].str();
                string node1 = match[3].str();
                string node2 = match[4].str();
                string number = match[5].str();
                string unit = match[6].str();
                double value;
                try {
                    value = stod(number);
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value\n";
                    continue;
                }
                string error = controller.handleError(name, circuit, value);
                if (error != "") {
                    cout << error;
                    continue;
                }
                if (!unit.empty()) {
                    switch (unit[0]) {
                        case 'G':
                            value *= 1e9;
                            break;
                        case 'M':
                            value *= 1e6;
                            break;
                        case 'k':
                        case 'K':
                            value *= 1e3;
                            break;
                        case 'm':
                            value *= 1e-3;
                            break;
                        case 'u':
                            value *= 1e-6;
                            break;
                        case 'n':
                            value *= 1e-9;
                            break;
                    }
                }
                cout << controller.addNewElement(node1, node2, name, value, circuit);
            } else if (regex_match(input, match, remove_element)) {
                if (match[1] != "R" && match[1] != "L" && match[1] != "C" && match[1] != "D") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                cout << controller.removeElement(match[1].str() + match[2].str(), circuit);
            } else if (regex_match(input, match, add_diode)) {
                string name = match[1].str() + match[2].str();
                string node1 = match[3];
                string node2 = match[4];
                string model = match[5];
                string error = controller.handleError(name, circuit, 0);
                if (error != "") {
                    cout << error;
                    continue;
                }
                if (match[1] != "D") {
                    cout << "Error: Element " << name << " not found in library\n";
                    continue;
                }
                if (model == "D" || model == "Z")
                    cout << controller.addDiode(node1, node2, name, model, circuit);
            } else if (regex_match(input, match, add_ground)) {
                if (match[1] != "GND") {
                    cout << "Error: Element " << match[1] << " not found in library\n";
                    continue;
                }
                cout << controller.addGround(match[2], circuit);
            } else if (regex_match(input, match, delete_ground)) {
                if (match[1] != "GND") {
                    cout << "Error: Element " << match[1] << " not found in library\n";
                    continue;
                }
                cout << controller.deleteGround(match[2], circuit);
            } else if (regex_match(input, match, show_nodes)) {
                controller.showNodes(circuit);
            } else if (regex_match(input, match, list)) {
                controller.list(circuit);
            } else if (regex_match(input, match, list_element)) {
                if (match[1] != "R" && match[1] != "C" && match[1] != "L" && match[1] != "D" && match[1] != "Z") {
                    cout << "Error: Element " << match[1] << " not found in library\n";
                    continue;
                }
                controller.listElement(match[1], circuit);
            } else if (regex_match(input, match, rename_node)) {
                string old_name = match[1];
                string new_name = match[2];
                string err;
                circuit->renameNode(old_name, new_name, err);
                cout << err;
            } else if (regex_match(input, match, show_details)) {
                controller.showCircuitDetails(circuit);
            } else if (regex_match(input, match, new_file)) {
                string address = match[1];
                ifstream fin(address, ios::in);

                if (!fin) {
                    cerr << "Error opening file!" << endl;
                    return;
                }

                cout << "lets start" << endl;

                string line;
                regex word_regex("\\S+");
                while (getline(fin, line)) {
                    auto words_begin = sregex_iterator(line.begin(), line.end(), word_regex);
                    auto words_end = sregex_iterator();

                    vector<string> words;

                    for (auto it = words_begin; it != words_end; ++it) {
                        words.push_back(it->str());
                    }

                    string name;
                    string node1;
                    string node2;
                    string tmpValue;
                    string tmpAmplitude;
                    string tmpFrequency;
                    string tmpOffset;
                    string tmpVi;
                    string tmpVf;
                    string tmpTD;
                    string tmpTR;
                    string tmpTF;
                    string tmpTOn;
                    string tmpPeriod;


                    string type = words[0];
                    if (type == "R" || type == "L" || type == "C" || type == "I" || type == "V" || type == "D" ||
                        type == "Z") {
                        name = words[1];
                        node1 = words[2];
                        node2 = words[3];
                        tmpValue = words[4];
                        if (node1 == "GND") {
                            controller.addGround(node1, circuit);
                        }
                        if (node2 == "GND") {
                            controller.addGround(node2, circuit);
                        }

                        double value = 0;

                        if (tmpValue[0] == '-') {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpValue.back() == 'G') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                            } else if (tmpValue.back() == 'M') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                            } else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                            } else if (tmpValue.back() == 'u') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                            } else if (tmpValue.back() == 'n') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                            } else if (tmpValue.back() == 'm') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                            } else {
                                value = stod(tmpValue);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }
                        if (type == "R" || type == "C" || type == "L") {
                            cout << controller.addNewElement(node1, node2, name, value, circuit);

                        } else if (type == "D" || type == "Z") {
                            string model;
                            if (type == "D") {
                                model = "D";
                            }
                            if (type == "Z") {
                                model = "Z";
                            }
                            cout << controller.addDiode(node1, node2, name, model, circuit);
                        } else if (type == "V") {
                            cout << controller.addDCVoltageSource(name, node1, node2, value, circuit);
                        } else if (type == "I") {
                            cout << controller.addDCCurrentSource(name, node1, node2, value, circuit);
                        }
                    }
                    if (type == "VSIN" || type == "ISIN") {
                        name = words[1];
                        node1 = words[2];
                        node2 = words[3];
                        tmpAmplitude = words[4];
                        tmpFrequency = words[5];
                        tmpOffset = words[6];
                        double amplitude;
                        double frequency;
                        double offset;
                        try {
                            if (tmpAmplitude.back() == 'G') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e9;
                            } else if (tmpAmplitude.back() == 'M') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e6;
                            } else if (tmpAmplitude.back() == 'k' || tmpAmplitude.back() == 'K') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e3;
                            } else if (tmpAmplitude.back() == 'u') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e-6;
                            } else if (tmpAmplitude.back() == 'n') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e-9;
                            } else if (tmpAmplitude.back() == 'm') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e-3;
                            } else {
                                amplitude = stod(tmpAmplitude);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpFrequency.back() == 'G') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e9;
                            } else if (tmpFrequency.back() == 'M') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e6;
                            } else if (tmpFrequency.back() == 'k' || tmpFrequency.back() == 'K') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e3;
                            } else if (tmpFrequency.back() == 'u') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e-6;
                            } else if (tmpFrequency.back() == 'n') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e-9;
                            } else if (tmpFrequency.back() == 'm') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e-3;
                            } else {
                                frequency = stod(tmpFrequency);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }
                        try {
                            if (tmpOffset.back() == 'G') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e9;
                            } else if (tmpOffset.back() == 'M') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e6;
                            } else if (tmpOffset.back() == 'k' || tmpOffset.back() == 'K') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e3;
                            } else if (tmpOffset.back() == 'u') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e-6;
                            } else if (tmpOffset.back() == 'n') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e-9;
                            } else if (tmpOffset.back() == 'm') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e-3;
                            } else {
                                offset = stod(tmpOffset);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        if (type == "VSIN") {
                            cout << controller.addSinusoidalVoltageSource(name, node1, node2, amplitude, frequency,
                                                                          offset, circuit);
                        }
                        if (type == "ISIN") {
                            cout << controller.addSinusoidalCurrentSource(name, node1, node2, amplitude, frequency,
                                                                          offset, circuit);
                        }
                    }

                    if (type == "VPULSE" || type == "IPULSE") {
                        name = words[1];
                        node1 = words[2];
                        node2 = words[3];
                        tmpVi = words[4];
                        tmpVf = words[5];
                        tmpTD = words[6];
                        tmpTR = words[7];
                        tmpTF = words[8];
                        tmpTOn = words[9];
                        tmpPeriod = words[10];
                        double Vi, Vf, TD, TR, TF, TOn, period;
                        try {
                            if (tmpVi.back() == 'G') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e9;
                            } else if (tmpVi.back() == 'M') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e6;
                            } else if (tmpVi.back() == 'k' || tmpVi.back() == 'K') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e3;
                            } else if (tmpVi.back() == 'u') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e-6;
                            } else if (tmpVi.back() == 'n') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e-9;
                            } else if (tmpVi.back() == 'm') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e-3;
                            } else {
                                Vi = stod(tmpVi);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpVf.back() == 'G') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e9;
                            } else if (tmpVf.back() == 'M') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e6;
                            } else if (tmpVf.back() == 'k' || tmpVf.back() == 'K') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e3;
                            } else if (tmpVf.back() == 'u') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e-6;
                            } else if (tmpVf.back() == 'n') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e-9;
                            } else if (tmpVf.back() == 'm') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e-3;
                            } else {
                                Vf = stod(tmpVf);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpTD.back() == 'G') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e9;
                            } else if (tmpTD.back() == 'M') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e6;
                            } else if (tmpTD.back() == 'k' || tmpTD.back() == 'K') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e3;
                            } else if (tmpTD.back() == 'u') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e-6;
                            } else if (tmpTD.back() == 'n') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e-9;
                            } else if (tmpTD.back() == 'm') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e-3;
                            } else {
                                TD = stod(tmpTD);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpTR.back() == 'G') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e9;
                            } else if (tmpTR.back() == 'M') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e6;
                            } else if (tmpTR.back() == 'k' || tmpTR.back() == 'K') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e3;
                            } else if (tmpTR.back() == 'u') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e-6;
                            } else if (tmpTR.back() == 'n') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e-9;
                            } else if (tmpTR.back() == 'm') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e-3;
                            } else {
                                TR = stod(tmpTR);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpTF.back() == 'G') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e9;
                            } else if (tmpTF.back() == 'M') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e6;
                            } else if (tmpTF.back() == 'k' || tmpTF.back() == 'K') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e3;
                            } else if (tmpTF.back() == 'u') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e-6;
                            } else if (tmpTF.back() == 'n') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e-9;
                            } else if (tmpTF.back() == 'm') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e-3;
                            } else {
                                TF = stod(tmpTF);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpTOn.back() == 'G') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e9;
                            } else if (tmpTOn.back() == 'M') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e6;
                            } else if (tmpTOn.back() == 'k' || tmpTOn.back() == 'K') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e3;
                            } else if (tmpTOn.back() == 'u') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e-6;
                            } else if (tmpTOn.back() == 'n') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e-9;
                            } else if (tmpTOn.back() == 'm') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e-3;
                            } else {
                                TOn = stod(tmpTOn);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpPeriod.back() == 'G') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e9;
                            } else if (tmpPeriod.back() == 'M') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e6;
                            } else if (tmpPeriod.back() == 'k' || tmpPeriod.back() == 'K') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e3;
                            } else if (tmpPeriod.back() == 'u') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e-6;
                            } else if (tmpPeriod.back() == 'n') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e-9;
                            } else if (tmpPeriod.back() == 'm') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e-3;
                            } else {
                                period = stod(tmpPeriod);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }
                        if (type == "VPULSE") {
                            cout << controller.addPulseVoltageSource(name, node1, node2, Vi, Vf, TD, TR, TF, TOn,
                                                                     period, circuit);
                        }
                        if (type == "IPULSE") {
                            cout << controller.addPulseCurrentSource(name, node1, node2, Vi, Vf, TD, TR, TF, TOn,
                                                                     period, circuit);
                        }
                    }
                    if (type == "E" || type == "G") {
                        name = words[1];
                        node1 = words[2];
                        node2 = words[3];
                        string controlNode1 = words[4];
                        string controlNode2 = words[5];
                        string tmpGain = words[5];
                        double gain;
                        Node *cn1 = circuit->getCreateNode(controlNode1);
                        Node *cn2 = circuit->getCreateNode(controlNode2);
                        bool cn1Connected = false;
                        bool cn2Connected = false;
                        for (auto element: circuit->getElements()) {
                            if (element->getFirstNode() == cn1 || element->getSecondNode() == cn1) {
                                cn1Connected = true;
                            }
                            if (element->getFirstNode() == cn2 || element->getSecondNode() == cn2) {
                                cn2Connected = true;
                            }
                        }

                        if (!cn1Connected || !cn2Connected) {
                            cout << "Error: " << controlNode1 << " and " << controlNode2 << " are not connected"
                                 << endl;
                            continue;
                        }
                        try {
                            if (tmpGain.back() == 'G') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e9;
                            } else if (tmpGain.back() == 'M') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e6;
                            } else if (tmpGain.back() == 'k' || tmpGain.back() == 'K') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e3;
                            } else if (tmpGain.back() == 'u') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-6;
                            } else if (tmpGain.back() == 'n') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-9;
                            } else if (tmpGain.back() == 'm') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-3;
                            } else {
                                gain = stod(tmpGain);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }
                        if (gain == 0) {
                            cout << "Error: Gain cannot be zero" << endl;
                            continue;
                        }
                        if (type == "E") {
                            cout << controller.addVCVS(name, node1, node2, controlNode1, controlNode2, gain, circuit);
                        }
                        if (type == "G") {
                            cout << controller.addVCCS(name, node1, node2, controlNode1, controlNode2, gain, circuit);
                        }
                    }

                    if (type == "H" || type == "F") {
                        name = words[1];
                        node1 = words[2];
                        node2 = words[3];
                        string cElement = words[4];
                        string tmpGain = words[5];
                        double gain;

                        Element *controlElement = nullptr;
                        for (auto &element: circuit->getElements()) {
                            if (element->getName() == cElement) {
                                controlElement = element;
                                break;
                            }
                        }
                        if (controlElement == nullptr) {
                            cout << "Error: " << cElement << " does not exist in the circuit" << endl;
                            continue;
                        }
                        try {
                            if (tmpGain.back() == 'G') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e9;
                            } else if (tmpGain.back() == 'M') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e6;
                            } else if (tmpGain.back() == 'k' || tmpGain.back() == 'K') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e3;
                            } else if (tmpGain.back() == 'u') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-6;
                            } else if (tmpGain.back() == 'n') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-9;
                            } else if (tmpGain.back() == 'm') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-3;
                            } else {
                                gain = stod(tmpGain);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }
                        if (gain == 0) {
                            cout << "Error: Gain cannot be zero" << endl;
                            continue;
                        }
                        if (type == "H") {
                            cout << controller.addCCVS(name, node1, node2, cElement, gain, circuit);
                        }
                        if (type == "F") {
                            cout << controller.addCCCS(name, node1, node2, cElement, gain, circuit);
                        }

                    }

                }

                cout << "reading file ended :)" << endl;
                fin.close();
            } else if (regex_match(input, match, addVCVS)) {
                string name = 'E'+match[1].str();
                string node1 = match[2].str();
                string node2 = match[3].str();
                string controlNode1 = match[4].str();
                string controlNode2 = match[5].str();
                string tmpValue = match[6].str();
                double gain;
                try {
                    if (tmpValue.back() == 'G') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                    } else if (tmpValue.back() == 'M') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                    } else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                    } else if (tmpValue.back() == 'u') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                    } else if (tmpValue.back() == 'n') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                    } else if (tmpValue.back() == 'm') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                    } else {
                        gain = stod(tmpValue);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (gain == 0) {
                    cout << "Error: Gain cannot be zero" << endl;
                    continue;
                }
                Node *cn1 = circuit->getCreateNode(controlNode1);
                Node *cn2 = circuit->getCreateNode(controlNode2);
                bool cn1Connected = false;
                bool cn2Connected = false;
                for (auto element: circuit->getElements()) {
                    if (element->getFirstNode() == cn1 || element->getSecondNode() == cn1) {
                        cn1Connected = true;
                    }
                    if (element->getFirstNode() == cn2 || element->getSecondNode() == cn2) {
                        cn2Connected = true;
                    }
                }

                if (!cn1Connected || !cn2Connected) {
                    cout << "Error: " << controlNode1 << " and " << controlNode2 << " are not connected";
                    continue;
                }
                cout << controller.addVCVS(name, node1, node2, controlNode1, controlNode2, gain, circuit);
            } else if (regex_match(input, match, addVCCS)) {
                string name = 'G'+match[1].str();
                string node1 = match[2].str();
                string node2 = match[3].str();
                string controlNode1 = match[4].str();
                string controlNode2 = match[5].str();
                string tmpValue = match[6].str();
                double gain;
                try {
                    if (tmpValue.back() == 'G') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                    } else if (tmpValue.back() == 'M') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                    } else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                    } else if (tmpValue.back() == 'u') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                    } else if (tmpValue.back() == 'n') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                    } else if (tmpValue.back() == 'm') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                    } else {
                        gain = stod(tmpValue);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (gain == 0) {
                    cout << "Error: Gain cannot be zero" << endl;
                    continue;
                }
                Node *cn1 = circuit->getCreateNode(controlNode1);
                Node *cn2 = circuit->getCreateNode(controlNode2);
                bool cn1Connected = false;
                bool cn2Connected = false;
                for (auto element: circuit->getElements()) {
                    if (element->getFirstNode() == cn1 || element->getSecondNode() == cn1) {
                        cn1Connected = true;
                    }
                    if (element->getFirstNode() == cn2 || element->getSecondNode() == cn2) {
                        cn2Connected = true;
                    }
                }

                if (!cn1Connected || !cn2Connected) {
                    cout << "Error: " << controlNode1 << " and " << controlNode2 << " are not connected";
                    continue;
                }
                cout << controller.addVCCS(name, node1, node2, controlNode1, controlNode2, gain, circuit);
            } else if (regex_match(input, match, addCCVS)) {
                string name = 'H'+match[1].str();
                string node1 = match[2].str();
                string node2 = match[3].str();
                string cElement = match[4].str();
                string tmpValue = match[5].str();
                double gain;
                try {
                    if (tmpValue.back() == 'G') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                    } else if (tmpValue.back() == 'M') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                    } else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                    } else if (tmpValue.back() == 'u') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                    } else if (tmpValue.back() == 'n') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                    } else if (tmpValue.back() == 'm') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                    } else {
                        gain = stod(tmpValue);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (gain == 0) {
                    cout << "Error: Gain cannot be zero" << endl;
                    continue;
                }

                Element *controlElement = nullptr;
                for (auto &element: circuit->getElements()) {
                    if (element->getName() == cElement) {
                        controlElement = element;
                        break;
                    }
                }
                if (controlElement == nullptr) {
                    cout << "Error: " << cElement << " does not exist in the circuit" << endl;
                    continue;
                }
                cout << controller.addCCVS(name, node1, node2, cElement, gain, circuit);
            } else if (regex_match(input, match, addCCCS)) {
                string name = 'F'+match[1].str();
                string node1 = match[2].str();
                string node2 = match[3].str();
                string cElement = match[4].str();
                string tmpValue = match[5].str();
                double gain;
                try {
                    if (tmpValue.back() == 'G') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                    } else if (tmpValue.back() == 'M') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                    } else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                    } else if (tmpValue.back() == 'u') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                    } else if (tmpValue.back() == 'n') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                    } else if (tmpValue.back() == 'm') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                    } else {
                        gain = stod(tmpValue);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (gain == 0) {
                    cout << "Error: Gain cannot be zero" << endl;
                    continue;
                }

                Element *controlElement = nullptr;
                for (auto &element: circuit->getElements()) {
                    if (element->getName() == cElement) {
                        controlElement = element;
                        break;
                    }
                }
                if (controlElement == nullptr) {
                    cout << "Error: " << cElement << " does not exist in the circuit" << endl;
                    continue;
                }
                cout << controller.addCCCS(name, node1, node2, cElement, gain, circuit);
            } else if (regex_match(input, match, preanalysischeck)) {
                controller.preAnalysisErrs(circuit);
            } else if (regex_match(input, match, exit)) {
                cout << "Bye Bye!\n";
                return;
            } else
                cout << "Syntax error\n";
        }
    }
};

int main(){
    View view;
    view.run();
    return 0;
}