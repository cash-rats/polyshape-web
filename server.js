import express from 'express';
import path from 'path';
const app = express();

app.use(express.static('build/client'));

app.get('*', (req, res) => {
  res.sendFile(path.join(path.__dirname, 'build/client/index.html'));
});

const PORT = process.env.PORT || 3000;
app.listen(PORT, () => {
  console.log(`Server running on port ${PORT}`);
});
