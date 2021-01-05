#include "types.h"

// useless function
PyObject *Match_new(PyTypeObject *type, PyObject *args, PyObject *kwargs){
    
    MATCH *self;
    
    self = (MATCH *)type->tp_alloc(type, 0);

    if (!self) return NULL;

    self->query[0] = '\0';
    self->subject[0] = '\0';
    self->cigar[0] = '\0';
    self->start = 0;
    self->end = 0;
    self->mq = 0;
    self->read_length = 0;
    self->percent_id = 0.0;
    self->paired = true;
    self->parity= true;
    self->mapped= false;
    self->orphan= true;
    self->multi= true;
    self->chimeric= true;
    self->singleton= true;
    self->w=0.0;

    return (PyObject *)self;
}


static int Match_init(MATCH *self, PyObject *args, PyObject *kwargs){
    static char* kwlist[] = {"start", "end" , "weight", "query", "cigar", "subject", "read_length", "percent_id", NULL};

    if(! PyArg_ParseTupleAndKeywords(args, kwargs, "iifsssif", kwlist,
           &self->start, &self->end, &self->w, &self->query, &self->cigar, &self->subject, &self->read_length, &self->percent_id))
        return -1;
    return 0;

}

static void Match_dealloc(MATCH *self){
    free(self->cigar);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyMemberDef Match_members[] = {
    {"start", T_UINT, offsetof(MATCH, start),0, "Match attribute"}, //T_UINT, T_INT
    {"end", T_UINT , offsetof(MATCH, end), 0, "Match attribute"},
    {"weight", T_FLOAT, offsetof(MATCH, w), 0, "Match attribute"},
    {"query", T_STRING , offsetof(MATCH, query), 0, "Match attribute"}, //string type are read_only after passing to python
    {"cigar", T_STRING , offsetof(MATCH, cigar), 0, "Match attribute"},
    {"subject", T_STRING , offsetof(MATCH, subject), 0, "Match attribute"},
    {"read_length", T_FLOAT, offsetof(MATCH, read_length), 0, "Match attribute"},
    {"percent_id", T_FLOAT , offsetof(MATCH, percent_id), 0, "Match attribute"},
    {NULL}
};

PyObject *Match_repr(MATCH *self) {
	if (self->mapped) {
		return PyUnicode_FromFormat("<Match> %s (%d) mapped to %s",
		                            self->query, self->read_length, self->subject);
	} else {
		return PyUnicode_FromFormat("<Match> %s (%d) was not aligned",
		                            self->query, self->read_length);
	}
}

MATCH *Match_cnew(PyTypeObject *type){
    MATCH *self;
    self = (MATCH *)type->tp_alloc(type,0);

    if (self) return self;
    else return NULL;

}

PyTypeObject MatchType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_sam_module.Match",       /* tp_name */
    sizeof(MATCH),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Match_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    (reprfunc)Match_repr,      /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "Match objects",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    0,//Match_methods,         /* tp_methods */
    Match_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Match_init,      /* tp_init */
    PyType_GenericAlloc,       /* tp_alloc */
    Match_new,                 /* tp_new */
};

unsigned int decode_cigar(MATCH* self){
    unsigned int read_len = 0;
    unsigned int aln_len = 0 ;

    string consume_ref =  "MDN=X";
	string consume_query = "MIS=X";
	
	string buffer = "";

    char * c;
    for (c = self->cigar; *c != '\0'; c++ ){
		if(isdigit(*c)){
			buffer = buffer + *c;
		} else {
			if(consume_ref.find(*c) != string::npos){
				aln_len += atoi(buffer.c_str());
			}
            if (consume_query.find(*c) != string::npos){
				read_len += atoi(buffer.c_str());
			}
			buffer = "";
        }
    }

    if (read_len == 0)
        PyErr_SetString(PyExc_ValueError, "alignment length calculated from CIGAR was zero.");

    self->read_length = read_len;

    return aln_len;


}

void update_end_and_read_length(MATCH * self){
    if (self->subject == "UNMAPPED")
        return ;

    unsigned int aln_len = decode_cigar(self);
    self->end = self->start + aln_len;
}