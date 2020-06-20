#include "types.h"

// useless function
PyObject *Match_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    
    MATCH *self;
    
    self = (MATCH *)type->tp_alloc(type,0);
    if (self != NULL){
        self->query = "";
        self->subject = "";
        self->cigar = "";
        self->start = 0;
        self->end = 0;
        self->mq = 0;
        self->read_length = 0;
        self->percent_id = 0.0;
        self->paired = true;
        self->parity= true;
        self->mapped= true;
        self->orphan= true;
        self->multi= true;
        self->chimeric= true;
        self->singleton= true;
        self->w=0.0;
    }
    
    return (PyObject *)self;
    
}

//useless function; initializer for Python; don't think this object will be instatiated in Python
//add this if needed
static int Match_init(MATCH *self, PyObject *args, PyObject *kwds){
    /*
    static char *kwlist[] = {"n", "m", NULL};
    
    if( ! PyArg_ParseTupleAndKeywords(args, kwds, "ii", kwlist, &self->n, &self->m))
       return -1;
    
    */
    static char *kwlist[] = {"start", "end" , "weight", "ref", "query", "cigar","subject" ,"read_length" , "percent_id", NULL};
    
    if(! PyArg_ParseTupleAndKeywords(args, kwds, "iifssssff",kwlist, &self->start, &self->end,&self->w,&self->ref,&self->query,&self->cigar,&self->subject,&self->read_length,&self->percent_id))
        return -1;
    return 0;
       
}

static void Match_dealloc(MATCH *self){
    Py_TYPE(self)->tp_free((PyObject*)self);
}

//static PyGetSetDef Car_getseters[]

static PyMemberDef Match_members[] = {
    {"start", T_UINT, offsetof(MATCH, start),0, "Match attribute"},//T_UINT, T_INT
    {"end", T_UINT , offsetof(MATCH, end), 0,"Match attribute"},
    {"weight", T_FLOAT, offsetof(MATCH, w), 0,"Match attribute"},
    {"ref", T_STRING, offsetof(MATCH, ref), 0,"Match attribute"},
    {"query", T_STRING , offsetof(MATCH, query), 0,"Match attribute"}, //string type are read_only after passing to python
    {"cigar", T_STRING , offsetof(MATCH, cigar), 0,"Match attribute"},
    {"subject", T_STRING , offsetof(MATCH, subject), 0,"Match attribute"},
    {"read_length", T_FLOAT, offsetof(MATCH, read_length), 0,"Match attribute"},
    {"percent_id", T_FLOAT , offsetof(MATCH, percent_id), 0,"Match attribute"},
    {NULL}
};

MATCH *Match_cnew(PyTypeObject *type){
    MATCH *self;
    self = (MATCH *)type->tp_alloc(type,0);

    if(self )

    return self;
}

PyTypeObject MatchType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_sam_module.Match",             /* tp_name */
    sizeof(MATCH),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Match_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,   /* tp_flags */
    "Match objects",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    0,//Match_methods,             /* tp_methods */
    Match_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Match_init,      /* tp_init */
    0,                         /* tp_alloc */
    Match_new,                 /* tp_new */
};