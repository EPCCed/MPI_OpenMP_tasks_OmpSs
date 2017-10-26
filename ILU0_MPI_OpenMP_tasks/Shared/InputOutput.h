#ifndef InputOutput

#define InputOutput 1

extern void ClearScreen ();

extern char ReadCharacter ();

extern char ReadAlfanum ();

extern unsigned int ReadUnsignedInt ();

extern int ReadInt ();

extern float ReadUnsignedReal ();

extern float ReadReal ();

extern double ReadUnsignedDouble ();

extern double ReadDouble ();

extern void WaitChar ();

extern int YesOrNot (char *str);

extern int SelectOption (int min, int max, char *str);

#define DoChanges 1

extern void ChangeBinary (void *c, int size);

extern void ReadBinaryFile (FILE *f, void *c, int size);

extern void WriteBinaryFile (FILE *f, void *c, int size);

extern FILE *OpenFile (char *name, char *attr);

extern void ReadStringFile (FILE *file, char *string, int length);

extern void PrintSpaceFile (FILE *file, int num);

void PrintTrace (void);

#endif

