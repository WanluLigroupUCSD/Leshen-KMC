---
name: User rejects Explore agent tool calls
description: User has rejected Explore subagent calls multiple times - may prefer direct analysis or different approach
type: feedback
---

User rejected Explore subagent calls twice when analyzing Zacros and Leshen-KMC source code.

**Why:** Unclear - possibly the agent tool is too slow, too verbose, or the user prefers a different workflow. May also be related to interrupting to save state first.
**How to apply:** Next time, consider reading files directly with Read tool instead of spawning Explore agents, or ask the user their preferred approach before launching heavy analysis agents.
